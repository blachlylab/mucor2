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
import time
import operator as op
import gzip
# Local imports
from annfield import annfield
class VCFwork:
    def __init__(self, filename,by_anno,merge_anno):
        self.fn=filename
        self.by_anno=by_anno
        self.merge_anno=merge_anno


def flatten(x):
    if type(x)==type([]):
        x=",".join(map(str,x))
    return x

def VCF_to_lists(vcfwork:VCFwork)->list:
    fn=vcfwork.fn
    combine=vcfwork.merge_anno
    by_anno=vcfwork.by_anno
    #start return array and import vcf
    ret=[]
    reader = vcf.VCFReader(filename=fn)
    #get format and info keys for header
    formats = list(reader.formats.keys())
    infos = sorted(list(reader.infos.keys()))
    anntrue=False
    #remove ANN field
    for x in ['ANN']:
        if x in infos:
            anntrue=True
            del infos[infos.index(x)]
    #Make header
    ret.append(["SAMPLE"] + ['CHROM', 'POS', 'REF', 'ALT', 'ID','FILTER']+\
               ['info_' + x for x in infos]+formats)
    if anntrue:
        ret[0]=ret[0]+list(next(annfield.decode("|||||||||||||||")).keys())
    #loop over every row in vcf
    for record in reader:
        info_row=[]
        for x in infos:
            info_row.append(str(record.INFO.get(x, "")))
        # if a row in the return array for every annotation
        # pass ANN to annfield
        anns=[]
        if anntrue:
            for x in annfield.decode(",".join(record.INFO.get("ANN",None))):
                anns.append(x)
        # for each sample make row to add
        fixed=[]
        for sample in record.samples:
            row = [sample.sample]
            # Identify correct alleles
            gentype=getattr(sample.data, 'GT', None).split("/")
            if gentype[0]=='.':
                continue
            for i,x in enumerate(gentype):
                gentype[i]=int(x)
            alts=[]
            for allele in gentype:
                if allele != 0:
                    alts.append(allele)
            #if row per annotation match annotation to allele
            for alt in alts:
                ann_list=anns
                if by_anno:
                    ann_list=[]
                    for x in anns:
                        if x["allele"]==str(record.ALT[alt-1]):
                            ann_list.append(x)
                    #combine annotations
                    ann_list=merge_annotations(ann_list,"feature_id")
                elif anntrue:
                    ann_list=merge_annotations(ann_list,"allele")
                #print(ann_list)
                row.append(record.CHROM)
                row.append(record.POS)
                row.append(record.REF)
                row.append(str(record.ALT[alt-1]))
                row.append(record.ID)
                record_filter=record.FILTER
                if record_filter == []: record_filter = ['PASS']
                # Format fields not present will simply end up "blank"
                # in the output
                row.append(flatten(record_filter or '.'))
                for x in info_row:
                    row.append(flatten(x))
                for x in formats:
                    row.append(flatten(getattr(sample.data, x, None)))

                #add either ANN column or a row per annotation
                #then submit to return array
                if anns==[]:
                    ret.append(row)
                else:
                    for ann in ann_list:
                        ret.append(row+ann)
    return ret


#merge annotations with identical feature_ids/transcripts
def merge_annotations(anns: list,merge_on:str)->list:
    dictkeys=list(anns[0].keys())
    keys={}
    #collect all unique transcripts
    for d in anns:
        if d[merge_on] in keys:
            for i,k in enumerate(dictkeys):
                keys[d[merge_on]][i].append(d[k])
        else:
            keys[d[merge_on]]=[]
            for i,k in enumerate(dictkeys):
                keys[d[merge_on]].append([])
                keys[d[merge_on]][i].append(d[k])
    #sys.exit(0)
    ret=[]
    for key in keys:
        anns=keys[key]
        if len(anns[0])>1:
            for i,val in enumerate(anns):
                anns[i]=";".join(np.unique(val))
        else:
            for i,key in enumerate(anns):
                anns[i]=anns[i][0]
        ret.append(anns)
    return ret

def alter_table(master:pd.DataFrame,config_args:dict):
    sub=master.copy()
    ad=sub["AD"].str.split(",",expand=True).astype(int)
    sub.insert(sub.columns.get_loc("AD"),"depth",ad.sum(axis=1))
    for id in ad.columns:
        id=int(id)
        col=""
        if id==0:
            col="AD_WT"
        else:
            col="AD_ALT"+str(id)
        sub.insert(sub.columns.get_loc("AD")+1+id,col,ad.iloc[:,id])
    qss=sub["QSS"].str.split(",",expand=True).astype(int)
    for id in qss.columns:
        id=int(id)
        col=""
        if id==0:
            col="QSS_WT"
        else:
            col="QSS_ALT"+str(id)
        sub.insert(sub.columns.get_loc("QSS")+1+id,col,qss.iloc[:,id])
    for id in ad.columns:
        id=int(id)
        col=""
        if id==0:
            col="QSS/AD_WT"
        else:
            col="QSS/AD_ALT"+str(id)
        sub.insert(sub.columns.get_loc("QSS_ALT"+str(ad.columns[-1]))+1+id,col,qss.iloc[:,id]/ad.iloc[:,id])

    return sub

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
    for i,op in enumerate(config_args["SELECT"]["NUM"]["OP"]):
        col1=config_args["SELECT"]["NUM"]["COL"][i][0]
        col2=config_args["SELECT"]["NUM"]["COL"][i][1]
        if type(col2)==int:
            if op==">":
                sub=sub[sub[col1]>col2]
    for col in config_args["SELECT"]["STAT"]:
        c=sub[col]
        sub=sub[c>(c.mean()-(2*c.std()))]
        sub=sub[c<(c.mean()+(2*c.std()))]
    return sub

#pivot dataframe
def pivot(master:pd.DataFrame,config_args:dict):
    master[config_args["PIVINDEX"]] = master[config_args["PIVINDEX"]].fillna("")
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
    write_dataframe(args.output+"_subset.tsv",sub,False)
    prot=pivot(sub,config_args)
    print('printing aggregate')
    write_dataframe(args.output+"_protsum.tsv",prot,False)


#import vcf files into a dictionary either single-threaded or multithreaded
def import_files(args):
    #Merge individual file dataframes into one dataframe
    #if multithreaded set up thread pool and execute
    if args.threads:
        print('importing files')
        pool=mp.Pool(processes=args.threads)
        tasks=[VCFwork(fn,args.by_anno,not args.dontcombine) for fn in args.input]
        res=pool.map_async(VCF_to_lists, tasks)
        while (True):
            if (res.ready()): break
            remaining = res._number_left
            print("\r{0:.2f}% imported".format(1-((remaining)/len(args.input))),end="")
            time.sleep(0.5)
        print()
        res=res.get()
        output=[]
        for i,p in enumerate(res):
            if i!=0:
                p=p[1::]
            output+=p
        print('building dataframe')
        return pd.DataFrame(output[1::],columns=output[0])
    #else process files one at a time
    else:
        print('importing files')
        frames=[]
        for i,vcf_file in enumerate(args.input):
            if i==0:
                frames+=VCF_to_lists(VCFwork(vcf_file,args.by_anno,not args.dontcombine))
            else:
                frames+=VCF_to_lists(VCFwork(vcf_file,args.by_anno,not args.dontcombine))[1::]
            print("{0:.2f}%".format((i+1)/len(args.input)),end="")
            print(" imported",end="\r")
        print('building dataframe')
        return pd.DataFrame(frames[1::],columns=frames[0])

def write_dataframe(fn:str,data:pd.DataFrame,hdf):
    if hdf:
        store = pd.HDFStore(fn)
        store.put(fn,data,format="table",data_columns=True)
        store.close()
    else:
        data.to_csv(fn,sep="\t")



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
    parser.add_argument('-nm','--nomaster',action='store_true',help="don't write master dataframe out to disk")
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
    if config_args["DOMASTERIDX"]:
        master_df.set_index(config_args["MASTERIDX"], inplace=True)
        if not args.nomaster:
            print('printing master')
            write_dataframe(args.output + "_master.tsv",master_df,False)
        master_df.reset_index(inplace=True)
    else:
        if not args.nomaster:
            print('printing master')
            write_dataframe(args.output + "_master.tsv",master_df,False)


    # Option 2, aggregate on gene symbol from annotations

    # TODO SHould we generate any summary metrics (row wise or group wise) ?
    print('filtering data')
    sub=alter_table(master_df,config_args)
    sub=filter_table(sub,config_args)
    print('aggregating data')
    aggregate_data(sub,args,config_args)
    # TODO output summary JSON

    valediction()

if __name__ == "__main__":
    main()
