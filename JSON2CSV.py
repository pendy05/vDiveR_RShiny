#import libraries
import argparse, sys
from os import path
import json
import pandas as pd
import numpy as np

#command-line options
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest="input", help='Absolute path to the DiMA output file in JSON format.',nargs='+')
parser.add_argument('-p', '--protein', dest="protein", help='The name of the given protein.', default=['Unknown Protein'],nargs='+')
parser.add_argument('-hs', '--hostspecie', dest="host", help='The name of the corresponding host.', default=['Unknown Host'],nargs='+')
parser.add_argument('-o', '--output',dest="output", help='Absolute path to the output file.', default="result.csv")
parser.add_argument('-l', '--lowsupport',dest="lowSupport", help='Threshold for low support', default=30)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

#function to convert json file to csv file
def json2csv(infile,protein,host):
    #read input file in a json format
    with open(infile,'r') as f:
        data = json.loads(f.read())

    #normalize semi-structured JSON data into a flat table
    df_unnested = pd.json_normalize(data['results'], record_path =['variants'], meta =['position', 'entropy', 'support','distinct_variants_count', 'distinct_variants_incidence'])
    
    #create a dataframe for entropy, count, kmer_types.incidence, sequence and position    
    df = df_unnested[df_unnested['motif_long']=='Index']

    #label positions with multiple index
    df.insert(1, "multiIndex", df.duplicated(subset=['position'], keep=False))

    #label positions with low support
    df.insert(1, "lowSupport", np.where(df['support'] < int(args.lowSupport), True, False))

    #let user to label data with a protein name as an independent column    
    df.insert(1, "proteinName", str(protein))  

    #let user to label data with a host name as an independent column    
    df.insert(1, "host", str(host))  

    #drop one of the Multiindex 
    df = df.drop_duplicates(subset=['position']) 

    #create a second dataframe for motif incidences and position
    df2 = pd.pivot_table(df_unnested, index='position', columns='motif_long', values='incidence', aggfunc=np.sum, fill_value=0)
    df2.reset_index(inplace=True)
    df2['Total'] = (100 - df2['Index'])

    #merge two dataframes on "position" column and rename columns
    save = df.merge(df2, on='position', how='right').rename(columns={"Index": "index.incidence", "Major": "major.incidence","Minor": "minor.incidence",
                         "Unique": "unique.incidence","Total": "totalVariants.incidence", "sequence" : "indexPeptide", "distinct_variants_incidence": "distinct_variants_incidence"})

    #change the order of columns
    column_names = ["proteinName", "position", "support","lowSupport", "entropy","indexPeptide", 
    "index.incidence","major.incidence", "minor.incidence","unique.incidence", "totalVariants.incidence","distinct_variants_incidence", "multiIndex", "host"]
    save = save.reindex(columns=column_names)

    #save output file in a csv format
    if path.exists(args.output): # if the file is existing, save without header
        print("The output is appended to the file \"%s\".\n"%args.output)
        save.to_csv(args.output, mode='a', header=None,index = False, sep=",")
    else:
        save.to_csv(args.output, mode='a', index = False, sep=",")

    #notify user
    print(f'Results of \"{infile}\" is successfully saved at {args.output}')

if __name__ == "__main__":
    #check if input files present in current directories
    for file in args.input:
        if path.exists(file) == False:
            sys.exit("ERROR: Input file \"%s\" is not found in current directory! "%file)
            
    df = pd.DataFrame()
    #one input file with one protein and host name
    if (len(args.input)==len(args.protein) and len(args.protein)==len(args.host) and len(args.host)==1):
        json2csv(args.input[0], args.protein[0], args.host[0]) #conversion occurs
        
    #check if same amount of inputs are inserted for all the three arguments (input, protein, and host)
    elif (len(args.input)==len(args.protein) and len(args.protein)==len(args.host)):
        for infile, protein, host in zip(args.input,args.protein,args.host):#conversion occurs
            json2csv(infile,protein,host) 
            
    #lesser input file(s) are input compared to protein and/ host names
    elif (len(args.input) < len(args.protein) or len(args.input) < len(args.host)):
        sys.exit("ERROR: Number of input file is lesser than the protein and/or the host names!\nEach of the input file(s) should associated with one protein name and one host name. Please check your input.\n")
    
    #more input file(s) are input compared to protein and/ host names   
    elif (len(args.input) > len(args.protein) or len(args.input) > len(args.host)):
        print("NOTE:\n=====")
        print("Input filenames detected: %s"%(args.input))
        print("Protein names detected:   %s"%(args.protein))
        print("Host names detected:      %s\n"%(args.host))
        print("Would you like to proceed with the missing value(s) replaced with \"unknown\"?(y/n):")
        reply= input().lower()
        
        if (reply == "y" or reply == "yes"):
            input_len = len(args.input)
            if len(args.protein) == len(args.host): #if both the length of protein and host name elements are the same, we can loop through them tgt
                while len(args.protein) < input_len and len(args.host) < input_len:
                    args.protein.append("Unknown Protein")
                    args.host.append("Unknown Host")
                for infile, protein, host in zip(args.input,args.protein,args.host): #conversion occurs
                    json2csv(infile,protein,host)   
                
            else:#else, loop through the protein and host name elements separately
                while len(args.protein) < input_len :
                    args.protein.append("Unknown Protein")
                while len(args.host) < input_len:
                    args.host.append("Unknown Host")
                for infile, protein, host in zip(args.input,args.protein,args.host): #conversion occurs
                    json2csv(infile,protein,host)
