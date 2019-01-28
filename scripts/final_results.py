#!/usr/bin/env python

import argparse
import sys
import os
import commands
import numpy
import subprocess
import time


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Generates a final results file where the calculated Bayes factors are included as columns into the original data file.')
    parser.add_argument('-i','--BS_data_file_name',action='store',dest='INPUT_FILE',type=str,required=True,help='The input data file name.')
    parser.add_argument('-f','--output_folder',action='store',dest='OUTPUT_FOLDER',type=str,required=False,help='The output folder where the results files will be stored.')
    parser.add_argument('-o','--output_file',action='store',dest='OUTPUT_FILE',type=str,required=True,help='Name for the final output file.')
    parser.add_argument('-b','--BF_file',action='store',dest='BF_FILE',type=str,required=True,help='The BF file name.')
    parser.add_argument('-w','--window_index_file',action='store',dest='WINDOW_INDEX',type=str,required=True,help='The window index file.')

    options = parser.parse_args()


    window_index=numpy.loadtxt(options.WINDOW_INDEX)

    BF=numpy.loadtxt(options.BF_FILE)

    BF_sorted=BF[BF[:,1].argsort()]

    BF_vector_for_output=numpy.chararray(len(window_index),itemsize=100)

    for i in range(0,len(window_index)):

            if window_index[i]!=0:
                BF_vector_for_output[i]=BF_sorted[int(window_index[i])-1][0]
            else:
                BF_vector_for_output[i]="*"


    output_file_id=options.OUTPUT_FILE.split('.')[0]    


    currenttime=time.localtime()
    temp_output_file="TEMP_%s_%s_%s_%s_%s_%s_%s.txt"%(output_file_id,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
    temp_output_file_sed_output="TEMP_%s_sed_%s_%s_%s_%s_%s_%s.txt"%(output_file_id,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])

    BF_vector_for_output.tofile("%s/%s"%(options.OUTPUT_FOLDER,temp_output_file),sep="\n")

    subprocess.call("sed \"s/'//g\" %s/%s > %s/%s"%(options.OUTPUT_FOLDER,temp_output_file,options.OUTPUT_FOLDER,temp_output_file_sed_output),shell=True)

    subprocess.call("paste %s %s/%s > %s/%s"%(options.INPUT_FILE,options.OUTPUT_FOLDER,temp_output_file_sed_output,options.OUTPUT_FOLDER,options.OUTPUT_FILE),shell=True)

    subprocess.call("rm %s/%s"%(options.OUTPUT_FOLDER,temp_output_file_sed_output),shell=True)
    subprocess.call("rm %s/%s"%(options.OUTPUT_FOLDER,temp_output_file),shell=True)
