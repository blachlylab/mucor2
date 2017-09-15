import argparse
import os
import sys
import vcf
import json
import collections
import copy
import pandas as pd
import multiprocessing as mp
import numpy as np
# Local imports
from annfield import annfield

def flatten(x):
    if type(x)==type([]):
        x=",".join(map(str,x))
    return x

def VCF_to_lists(fn:str,combine:bool,by_anno:bool)->list:
    #start return array and import vcf
    ret=[]
    f=open(fn,'r')
    reader = vcf.VCFReader(f)
    #get format and info keys for header
    formats = list(reader.formats.keys())
    infos = sorted(list(reader.infos.keys()))
    #remove ANN field
    for x in ['ANN']:
        if x in infos:
            del infos[infos.index(x)]
    #Make header
    ret.append(["SAMPLE"] + ['CHROM', 'POS', 'REF', 'ALT', 'ID','FILTER']+\
               ['info.' + x for x in infos]+formats)
    #loop over every row in vcf
    for record in reader:
        info_row = [record.INFO.get(x, None) for x in infos]
        # if a row in the return array for every annotation
        # pass ANN to annfield
        anns=[x for x in annfield.decode(",".join(record.INFO.get("ANN",None)))]
        #append keys to header
        if len(ret)==1:
            ret[0]=ret[0]+list(anns[0].keys())
        # for each sample make row to add
        fixed=[]
        for sample in record.samples:
            row = [sample.sample]
            # Identify correct alleles
            gentype=getattr(sample.data, 'GT', None).split("/")
            if gentype[0]=='.':
                continue
            gentype=[int(x) for x in gentype]
            alts=[allele for allele in gentype if allele != 0]
            #if row per annotation match annotation to allele
            for alt in alts:
                ann_list=anns
                if by_anno:
                    ann_list=[x for x in anns if x["allele"]==str(record.ALT[alt-1])]
                    #if combining annotations with same feature_id
                    #combine annotations
                    if combine:
                        ann_list=merge_annotations(ann_list,"feature_id")
                else:
                    ann_list=merge_annotations(ann_list,"allele")
                fixed = [record.CHROM, record.POS, record.REF, str(record.ALT[alt-1]),
                         record.ID]
                record_filter=record.FILTER
                if record_filter == []: record_filter = ['PASS']
                # Format fields not present will simply end up "blank"
                # in the output
                row += fixed
                row += [flatten(record_filter or '.')]
                row += [flatten(x) for x in info_row]
                row += [flatten(getattr(sample.data, x, None)) for x in formats]

                #add either ANN column or a row per annotation
                #then submit to return array
                for ann in ann_list:
                    ret.append(row+[ann[key] for key in ann])
    f.close()
    return ret

#appends provided list with a dictionaries of vcf rows


#get metadata from vcf
def getmeta(vcf_reader:vcf.model):
    metadata = collections.OrderedDict()
    metadata['infos'] = \
        { k:vcf_reader.infos[k]._asdict() for k in vcf_reader.infos}
    metadata['formats'] = \
        { k:vcf_reader.formats[k]._asdict() for k in vcf_reader.formats}
    return metadata


#merge annotations with identical feature_ids/transcripts
def merge_annotations(anns: list,merge_on:str)->list:
    transcripts=[]
    #collect all unique transcripts
    for d in anns:
        transcripts.append(d[merge_on])
    if len(transcripts)==len(np.unique(transcripts)):
        return anns
    transcripts=np.unique(transcripts)
    ret=[]
    #for each transcript collect keys
    # of all annotations that contain the specified
    #transcript
    for ann in transcripts:
        ds=[]
        keys=[]
        for d in anns:
            keys=keys+list(d.keys())
            if d[merge_on]==ann:
                ds.append(d)
        keys=np.unique(keys)
        #merge all dictionaries with ann
        #merge based on type
        d={}
        for k in keys:
            k_list=[]
            for d in ds:
                if k in d:
                    k_list.append(d[k])
            k_list=np.unique(k_list)
            if type(k_list[0])==str:
                d[k]=";".join(k_list)
            elif type(k_list[0])==list or len(k_list)>1:
                d[k]=";".join(str(e)for e in k_list)
            else:
                d[k]=k_list[0]
        ret.append(d)
    return ret


#filters dataframe based on config file
def filter_table(master:pd.DataFrame,config_args:dict):
    sub=master.copy()
    select={}
    drop={}
    #Replace empty strings with NaN
    sub.replace([""], np.nan, inplace=True)
    #Get all selects
    for item in config_args["SELECT"]["EXACT"]:
        if item not in select:
            select[item]=[]
        d=config_args["SELECT"]["EXACT"][item]
        for val in d:
            select[item].append("^"+val+"$")
    for item in config_args["SELECT"]["CONTAINS"]:
        if item not in select:
            select[item]=[]
        select[item]+=config_args["SELECT"]["CONTAINS"][item]
    #Get all drops
    for item in config_args["DROP"]["EXACT"]:
        if item not in select:
            drop[item]=[]
        d=config_args["DROP"]["EXACT"][item]
        for val in d:
            drop[item].append("^"+val+"$")
    for item in config_args["DROP"]["CONTAINS"]:
        if item not in select:
            drop[item]=[]
        drop[item]+=config_args["DROP"]["CONTAINS"][item]
    #drop rows that are null for columns in select list and dropnull list
    for item in config_args["DROPNULL"]+list(select):
        sub=sub[sub[item].notnull()]
    print(drop)
    print(select)
    #for columns selected keep rows with specified values in columns
    for key,val in select.items():
        sub=sub[sub[key].str.match("|".join(val))]
    #for columns selected drop rows with specified values in columns
    for key,val in drop.items():
        sub=sub[~sub[key].str.match("|".join(val))]
    return sub


#pivot dataframe
def pivot(master:pd.DataFrame,config_args:dict):
    sub = master[config_args["PIVINDEX"]+config_args["PIVON"]+config_args["PIVVAL"]]
    if config_args["PIVFUNC"]=="":
        return pd.pivot_table(sub,index=config_args["PIVINDEX"],columns=config_args["PIVON"],
                        fill_value=0,values=config_args["PIVVAL"])
    else:
        return pd.pivot_table(sub,index=config_args["PIVINDEX"],columns=config_args["PIVON"],
                        fill_value=0,values=config_args["PIVVAL"],aggfunc=config_args["PIVFUNC"])


#confirm input files exist
def validate(args):
    if not os.path.exists(args.config):
        print("{} does not exist".format(args.config))
        sys.exit(1)
    # Prevalidate so as not to interrupt long-running procedure
    for vcf_file in args.input:
        # Validate file exists
        if not os.path.exists(vcf_file):
            print("{} does not exist".format(vcf_file))
            sys.exit(1)


#TODO: Make a list of pivot functions for use
def set_piv_func(args):
    func={"count":np.count_nonzero}


#prints pivoted dataframe
def aggregate_data(sub:pd.DataFrame,args:argparse.ArgumentParser,config_args:dict):
    sub.to_csv(args.output+"_subset.tsv",sep="\t")
    prot=pivot(sub,config_args)
    print('printing aggregate')
    prot.to_csv(args.output+"_protsum.tsv",sep="\t")


#import vcf files into a dictionary either single-threaded or multithreaded
def import_files(args):
    #Merge individual file dataframes into one dataframe
    #if multithreaded set up thread pool and execute
    if args.threads:
        print('importing files')
        pool=mp.Pool(processes=args.threads)
        res=[pool.apply_async(VCF_to_lists, args=(fn,not args.dontcombine,args.by_anno)) for fn in args.input]
        output=[]
        for i,p in enumerate(res):
            if i==0:
                p=p.get()
            else:
                p=p.get()[1::]
            print("{0:.2f}%".format((i+1)/len(args.input)),end="")
            print(" imported",end="\r")
            output+=p
        print('building dataframe')
        return pd.DataFrame(output[1::],columns=output[0])
    #else process files one at a time
    else:
        print('importing files')
        frames=[]
        for i,vcf_file in enumerate(args.input):
            if i==0:
                frames+=VCF_to_lists(vcf_file,not args.dontcombine,args.by_anno)
            else:
                frames+=VCF_to_lists(vcf_file,not args.dontcombine,args.by_anno)[1::]
            print("{0:.2f}%".format((i+1)/len(args.input)),end="")
            print(" imported",end="\r")
        print('building dataframe')
        return pd.DataFrame(frames[1::],columns=frames[0])


#TODO make config

def greeting():
    print("Mucor 3")


def valediction():
    print("Have a nice day.")


def main():
    greeting()
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', action='append', help="Annotation field to include in output (default: all)")
    parser.add_argument('-i', '--input',nargs="+", help="Input VCF file(s) (can specify more than one with wildcard)")
    parser.add_argument('-o', '--output', help="Output file base name (exclude extension)")
    parser.add_argument('-t', '--threads', help="Number of threads to generate", type=int)
    parser.add_argument('-b', '--by_anno', action='store_true',help='create new row for each annotation')
    parser.add_argument('-c', '--config', default="config.json",help="Specify json config file")
    parser.add_argument('-mc', '--makeconfig',action='store_true', help="forms example config file")
    parser.add_argument('-dc','--dontcombine',action='store_true',
                        help="don't combine annotations with identical feature ids")
    args = parser.parse_args()
    if args.makeconfig:
        make_config(args)
    validate(args)
    #import config settings
    config_args=json.load(open(args.config))
    master_df=import_files(args)
    # subset on columns included in --annotation
    if args.annotation:
        # TODO print diagnostic pre subset
        master_df = master_df[ args.annotation ]
        # TODO print diagnostic post subset
    # TODO pivot and/or aggregate
    # Option 1, aggregate on location (+- X>Y change)
    print('printing master')
    if config_args["DOMASTERIDX"]:
        master_df.set_index(config_args["MASTERIDX"], inplace=True)
        master_df.to_csv( args.output + "_master.tsv",sep="\t")
        master_df.reset_index(inplace=True)
    else:
        master_df.to_csv( args.output + "_master.tsv",sep="\t")


    # Option 2, aggregate on gene symbol from annotations

    # TODO SHould we generate any summary metrics (row wise or group wise) ?
    print('filtering data')
    sub=filter_table(master_df,config_args)
    print('aggregating data')
    aggregate_data(sub,args,config_args)
    # TODO output summary JSON

    valediction()

if __name__ == "__main__":
    main()
