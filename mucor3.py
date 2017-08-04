import argparse
import os
import sys
import vcf
import json
import collections
import copy

import pandas as pd
import multiprocessing as mp

# Local imports
from annfield import annfield

#parses dictionary from vcf file into a dataframe
def consume_vcf(fn: str, by_anno: bool):
    frames = []
    for record in vcf_transform(fn):
        #if -b flag is used a row is created for each annotation in each row of
        #the vcf file
        if by_anno:
            for i, ann in enumerate(record["INFO"]["ANN"]):
                frames.append(copy.deepcopy(record))
                id=len(frames)-1
                frames[id].update(ann)
                frames[id]["ann_id"]=i
                del(frames[id]["INFO"]["ANN"])
        #else one row in the vcf file is one row in the dataframe
        else:
            row=flatten_dict_to_str(record)
            frames.append(row)
    return pd.DataFrame(frames)

#flatten lists in the dictionary into strings
def flatten_dict_to_str(di: dict):
    for key in di.keys():
        if type(di[key])==list:
            di[key]=", ".join(str(e)for e in di[key])
    return di

#strips the extension and the parent directories froma file path
def strip_filename(fn: str):
    fn=fn[fn.rfind("/")+1::]
    return fn[0:fn.find(".")]


def vcf_transform(filename: str, metafile: str=None, json_array: bool=False) -> None:
    """Parse a VCF file into a metadata component (header) and individual JSON objects. Each row (and for multisample VCFs, sample within a row) is written as a JSON array.

    :param filename: filename of VCF to open
    :type filename: str.
    :param metafile: filename of header metadata to write
    :type metafile: str
    :returns: None
    :raises: Nothing

    Args:
        filename (str): filename of VCF to open
        metafile (str): filename of header metadata to write

    Returns:
        NoneType: Nothing!
    """

    with open(filename,'r') as fi:
        vcf_reader = vcf.Reader(fi)

        # Collect the header metadata from VCF file
        metadata = collections.OrderedDict()
        metadata['filename'] = vcf_reader.filename
        # General VCF metadata fields
        metadata['metadata'] = vcf_reader.metadata

        # Defined VCF header fields (not incl contig=)
        # vcf.Reader.Filter, .Format, .Info are not serializable
        # into JSON by default (they are namedtuples),
        # so must convert them to dict first
        metadata['filters'] = \
            { k:vcf_reader.filters[k]._asdict() for k in vcf_reader.filters}
        metadata['formats'] = \
            { k:vcf_reader.formats[k]._asdict() for k in vcf_reader.formats}
        metadata['infos'] = \
            { k:vcf_reader.infos[k]._asdict() for k in vcf_reader.infos}

        # Write the header metadata from VCF file (if requested)
        if metafile:
            with open(metafile, 'w') as fo:
                fo.write(json.dumps(metadata, indent=4))

        # Iterate through the rows of VCF and print individual variant records
        if json_array: print("[")
        for record in vcf_reader:
            yield get_vcf_record(record, metadata, json_array)
        if json_array: print("    {}\n]") # deal with terminal comma

def get_vcf_record(record: vcf.model._Record, metadata: dict, json_array=False) -> None:
    """Parse a PyVCF record, which may contain multiple samples, and emit JSON"""

    # Reset with every VCF record (row)
    record_dict = {}
    record_dict['type'] = 'variant_vcf'
    record_dict['CHROM'] = record.CHROM
    record_dict['POS'] = record.POS
    record_dict['ID'] = record.ID
    record_dict['REF'] = record.REF
    record_dict['ALT'] = str(record.ALT[0])
    record_dict['QUAL'] = record.QUAL
    record_dict['FILTER'] = record.FILTER

    # seq: URI
    # TODO: write function; be smarter
    if 'hg19' in metadata['metadata']['reference']:
        record_dict['sequri'] = "seq:hg19/" + record.CHROM + ":" + str(record.POS)

    # There is a bug in PyVCF 0.6.8:
    # "PASS" is reported as empty list (as if value were '.')
    if record.FILTER == []: record_dict['FILTER'] = ['PASS']

    # INFO field is common to vcf record;
    # processing should be in outer loop (now this function)
    record_dict['INFO'] = {k:v for k,v in record.INFO.items() if v is not None}
    if 'ANN' in record.INFO:
        anno_raw = record.INFO['ANN']
        record_dict['INFO']['ANN'] = [x for x in annfield.decode(','.join(anno_raw))]
    if 'LOF' in record.INFO:
        # TODO decode LOF field (snpEff)
        pass
    if 'NMD' in record.INFO:
        # TODO decode NMD field (snpEff)
        pass

    for s in record.samples:
        # Compile sample level data (terminal column(s) in VCF)

        # VCF spec (4.1, sec 1.4.2) says that if genotype fields are present,
        # then the same types of data must be present for all samples
        sample_dict = {}
        sample_dict['sample'] = s.sample
        for key in s.data._fields:
            sample_dict[key + "|" + metadata['formats'][key]['desc'] ] = s[key]

        # Now merge the VCF record (unchanged for every sample within a row)
        # with the sample-specific data. Separating record specific data
        # provides massive speed boost in multisample VCF processing
        # NB: ** dict unpacking operator req Python 3.5 or greater
        denormalized_dict = \
            {k:v for k,v in { **record_dict, **sample_dict }.items() if v is not None}

        if json_array:
            print("    " + json.dumps(denormalized_dict) + ",")
        else:
            #print(json.dumps(denormalized_dict))
            return denormalized_dict

        # for k in ext.data._fields:
        #     dictionary[k + "|" + vcf_file.formats[k].desc] = ext[k]
        # dictionary = {k:v for k,v in dictionary.items() if v is not None}
        # print (json.dumps(dictionary))
        # dictionary = {}


def greeting():
    print("Mucor 3")

def valediction():
    print("Have a nice day.")

def main():
    greeting()

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', action='append', help="Annotation field to include in output (default: all)")
    parser.add_argument('-i', '--input',nargs="+", help="Input VCF file (can specify more than once)")
    parser.add_argument('-o', '--output', help="Output file base name (exclude extension)")
    parser.add_argument('-t', '--threads', help="Number of threads to generate", type=int)
    parser.add_argument('-b', '--by_anno', action='store_true',help='create new row for each annotation')



    args = parser.parse_args()

    # Prevalidate so as not to interrupt long-running procedure
    for vcf_file in args.input:
        # Validate file exists
        if not os.path.exists(vcf_file):
            print("{} does not exist".format(vcf_file))
            sys.exit(1)
    master_df=pd.DataFrame()
    #Merge individual file dataframes into one dataframe
    #if multithreaded set up thread pool and execute
    if args.threads:
        pool=mp.Pool(processes=args.threads)
        res=[pool.apply_async(consume_vcf, args=(fn,args.by_anno,)) for fn in args.input]
        output=[p.get() for p in res]
        master_df = pd.concat(output)
    #else process files one at a time
    else:
        frames=[]
        for vcf_file in args.input:
            vcf_df = consume_vcf(vcf_file,args.by_anno)
            frames.append(vcf_df)
        master_df = pd.concat(frames)
    # subset on columns included in --annotation
    if args.annotation:
        # TODO print diagnostic pre subset
        master_df = master_df[ args.annotation ]
        # TODO print diagnostic post subset
    master_df.to_csv( args.output + "_master.tsv", sep="\t" )
    # TODO pivot and/or aggregate
    # Option 1, aggregate on location (+- X>Y change)
    grouped = master_df.groupby(by=['CHROM','POS','REF','ALT'], as_index=True, sort=True)

    # Option 2, aggregate on gene symbol from annotations

    # TODO SHould we generate any summary metrics (row wise or group wise) ?
    #grouped=grouped.aggregate()
    #grouped.to_csv(args.output +".tsv",sep="\t")

    # TODO output summary JSON

    valediction()

if __name__ == "__main__":
    main()
