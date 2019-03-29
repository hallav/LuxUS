import numpy
import scipy
import pickle
import argparse
import json



#Sigmoid function
def sigmoid(x):
    if x.shape==():
        out=1.0/(1.0+numpy.exp(-x))
    else:
        out = numpy.zeros([x.shape[0]])
        for i in range(0,x.shape[0]):
            out[i]=1.0/(1.0+numpy.exp(-x[i]))
    return out


def generate_coordinates(Wrange,cytosines,predictors):
    #Generate random coordinates
    coords=numpy.zeros([1,cytosines])	
    while numpy.unique(coords).shape!=(cytosines,):
        coords = numpy.random.randint(1,Wrange,cytosines)

    coords = numpy.sort(coords)
	

    #Replicate the coordinates to be of the right dimension
    rep_coords=numpy.zeros([cytosines*predictors])
    cind = 0 #Current index
    for coord in coords:
        for j in range(0,predictors):
            rep_coords[cind] = coord
            cind =cind+1
  
  
    return rep_coords


def generate_B(sigmaB2,differential,mean_B):
    #mean_B is a vector of length N_predictors

    #The number of predictors is limited to two for this function

    if differential:
        muB=[mean_B[0],mean_B[1]]
    else:
        muB=[mean_B[0],0]

    B=numpy.zeros(2)
    B[0] = numpy.random.normal(muB[0],numpy.sqrt(sigmaB2))
    B[1] = muB[1]

    return B


def generate_D(n_predictors,n_replicates):

    D=numpy.zeros((n_replicates*2,n_predictors))
    ind=0
    for i in range(0,n_replicates):
        for condition in [0,1]:
            if condition==1:
                D_row=[1,1] #changed on 28.3.2018 
            else: 
                D_row=[1,0]
            D[ind,:]=D_row
            ind+=1

    return D




def generate_indicator_vectors(n_cytosines,n_replicates):

    ind_cyt=numpy.zeros(n_replicates*2*n_cytosines)
    ind_rep=numpy.zeros(n_replicates*2*n_cytosines)

    for i in range(0,n_cytosines):
        ind_cyt[range(i*2*n_replicates,(i+1)*2*n_replicates)]=i
        ind_rep[range(i*2*n_replicates,(i+1)*2*n_replicates)]=range(0,2*n_replicates)

    return ind_cyt,ind_rep



def generate_u_R(N_replicates,sigmar2):

    u_R=numpy.random.normal(loc=0,scale=numpy.sqrt(sigmar2),size=N_replicates*2)

    return u_R


def generate_u_C(n_cytosines,sigmaC2,l,coords):

    #Form SigmaC covariance matrix
    SigmaC=numpy.zeros((n_cytosines,n_cytosines))
    for i in range(0,n_cytosines):
        for j in range(0,n_cytosines):
            if i==j:
                SigmaC[i,j]=sigmaC2
            else:
                SigmaC[i,j]=sigmaC2*numpy.exp(-abs(coords[i]-coords[j])/pow(l,2))

    u_C=numpy.random.multivariate_normal(numpy.zeros(n_cytosines),SigmaC)

    return u_C  


def generate_Y(X,B,Z_C,u_C,Z_R,u_R,sigmae2,n_replicates,n_cytosines):

    #Z_C and Z_R are indicator vectors, not actual design matrices

    #Generate Y
    Y=numpy.zeros(n_cytosines*n_replicates*2)

    for i in range(0,n_replicates*2*n_cytosines):
        Y[i]=numpy.random.normal(numpy.dot(X[i],B)+u_C[Z_C[i]]+u_R[Z_R[i]],sigmae2)

    return Y


def generate_counts(n_reads,Y,n_replicates,n_cytosines,SE,BSE,BSBE):

    theta=numpy.zeros(n_replicates*2*n_cytosines)
    counts=numpy.zeros(n_cytosines*2*n_replicates)
    for i in range(0,n_cytosines*2*n_replicates):
        theta[i]=sigmoid(Y[i])
        p = theta[i]*((1-SE[i])*(1-BSE[i])+SE[i]*BSE[i])+(1-theta[i])*((1-SE[i])*(1-BSBE[i])+SE[i]*BSBE[i])
        counts[i] = numpy.random.binomial(n_reads,p)

    return counts,theta





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulates bisulphite sequencing data from the LuxUS model.')
    parser.add_argument('-e','--sigmae2',action='store',dest='sigmae2',type=float,required=True,help='Variance for residual term e.')
    parser.add_argument('-q','--reads',action='store',dest='reads',type=int,required=True,help='Number of BS calls i.e. (total) reads.')
    parser.add_argument('-r','--replicates',action='store',dest='replicates',type=int,required=True,help='Number of replicates for cases and controls, the total number of samples will be 2x replicates.')
    parser.add_argument('-c','--cytosines',action='store',dest='N_cytosines',type=int,required=True,help='Number of cytosines.')
    parser.add_argument('-l','--lengthscale',action='store',dest='l',type=int,required=True,help='Lengthscale parameter for covariance matrix for the simulations.')
    parser.add_argument('-w','--width',action='store',dest='width',type=int,required=True,help='Width of the genomic region where the N_cytosines lie (in basepairs).')
    parser.add_argument('-n','--datasets',action='store',dest='N_datasets',type=int,required=True,help='Number of datasets to be created for both cases (no differential methylation, differential methylation). Total number of resulting datasets is 2xN_datasets.')
    parser.add_argument('-f','--folder',action='store',dest='folder',type=str,required=True,help='Folder where to store the generated data files. Format /level1/level2')
    parser.add_argument('-g','--ID',action='store',dest='ID',type=str,required=True,help='Identifier used in the resulting file names.')
    parser.add_argument('-v','--sigmaB2',action='store',dest='sigmaB2',type=float,required=True,help='Variance for B, which is inserted into covariance matrix diagonal (SigmaB). Used for simulations.')
    parser.add_argument('-d','--sigmar2',action='store',dest='sigmar2',type=float,required=True,help='sigma_r2 used for simulations.')
    parser.add_argument('-a','--sigmac2',action='store',dest='sigmac2',type=float,required=True,help='sigma_c2 used for simulations.')
    parser.add_argument('-m','--mean_B',action='store',dest='mean_B',type=str,required=True,help='Mean for B0 (coefficient for the intercept term) and value for B1 (coefficient for case/control covariate) for the cases with differential methylation. Should be given as string in format [B1, B2]')
    parser.add_argument('-j','--sigmaB2_stan',action='store',dest='sigmaB2_stan',type=float,required=True,help='Variance for B, the value which will be used in LuxUS analysis. (Can be different from what is used for simulations).')
    parser.add_argument('-x','--save_LuxUS',action='store',dest='save_LuxUS',type=int,required=True,help='0 or 1 indicating whether the seimulated data should be saved in a format supported by LuxUS.')
    parser.add_argument('-y','--save_proportion_table',action='store',dest='save_proportion_table',type=int,required=True,help='0 or 1 indicating whether the seimulated data should be saved in proportion table format, which can be used with eg. Radmeth. Also a design matrix is saved into a text file.')
    parser.add_argument('-z','--save_LuxUS_sep',action='store',dest='save_LuxUS_sep',type=int,required=True,help='0 or 1 indicating whether the seimulated data should be saved in a format supported by LuxUS analysis where each cytosine is analysed separately.')

    options = parser.parse_args()

    mean_B = numpy.array(json.loads(options.mean_B))

    N_reads=options.reads
    N_replicates=options.replicates
    N_cytosines=options.N_cytosines
    N_predictors=2 #Number of predictors is fixed in this script
    l=options.l
    width=options.width

    N_datasets=options.N_datasets
    outputFolder=options.folder
	
    #variances and standard deviations
    sigmaB2 = options.sigmaB2
    sigmar2 = options.sigmar2
    sigmac2 = options.sigmac2
    sigmae2 = options.sigmae2


    mu_B = numpy.zeros([N_cytosines*N_predictors])
    mu_B_sep = numpy.zeros([N_predictors])


    fileID="%s_C%s_Q%s_R%s_W%s"%(options.ID,N_cytosines,N_reads,N_replicates,width)

    #Prior parameters for variance terms (to be used in the analysis with LuxUS), not used in simulations!

    a_C=6
    b_C=3

    a_E=5
    b_E=5

    a_R=98
    b_R=143


    for i in range(0,N_datasets):
        print("Dataset %s"%(i))	

        #Simulating experimental parameters
        bsEff_true=numpy.random.beta(99,1,size=N_replicates*2)
        bsbEff_true=numpy.random.beta(1,999,size=N_replicates*2)
        seqErr_true=numpy.random.beta(1,999,size=N_replicates*2)

        bsEff_true_rep = numpy.kron(numpy.ones((N_cytosines)),bsEff_true)
        bsbEff_true_rep = numpy.kron(numpy.ones((N_cytosines)),bsbEff_true)
        seqErr_true_rep = numpy.kron(numpy.ones((N_cytosines)),seqErr_true)

        bsbEff_mean=1.0/(1.0+999.0)
        seqErr_mean=1.0/(1.0+999.0)

        simulate_conversion_efficiency=numpy.random.binomial(10,numpy.tile(bsEff_true,(3066,1)))
        bsEff_estimate=numpy.sum(simulate_conversion_efficiency,axis=0)/(10.0*3066)

        print("The simulated values for (bsEff,bsbEff,seqErr)=(%s,%s,%s)"%(bsEff_true,bsbEff_true,seqErr_true))

        bsEff_rep = numpy.kron(numpy.ones((N_cytosines)),bsEff_estimate)
        bsbEff_rep = numpy.kron(numpy.ones((N_cytosines*N_replicates*2)),bsbEff_mean)
        seqErr_rep = numpy.kron(numpy.ones((N_cytosines*N_replicates*2)),seqErr_mean)
        reads_rep = numpy.kron(numpy.ones((N_cytosines*N_replicates*2)),N_reads)

        bsEff_rep_sep = bsEff_estimate.copy()
        bsbEff_rep_sep = numpy.kron(numpy.ones((N_replicates*2)),bsbEff_mean)
        seqErr_rep_sep = numpy.kron(numpy.ones((N_replicates*2)),seqErr_mean)
        reads_rep_sep = numpy.kron(numpy.ones((N_replicates*2)),N_reads)

        #Generate coordinates
		
        coordinates = generate_coordinates(width,N_cytosines,N_predictors)

        for diff in [0,1]:

	
            print("Differential=%s"%(diff))

            B=generate_B(sigmaB2,diff,mean_B)
            D=generate_D(N_predictors,N_replicates)
            Z_C,Z_R = generate_indicator_vectors(N_cytosines,N_replicates)
            u_R = generate_u_R(N_replicates,sigmar2)
            u_C = generate_u_C(N_cytosines,sigmac2,l,coordinates[::N_predictors])
            Y = generate_Y(numpy.tile(D,(N_cytosines,1)),B,Z_C.astype(int),u_C,Z_R.astype(int),u_R,sigmae2,N_replicates,N_cytosines)
            counts,thetas = generate_counts(N_reads,Y,N_replicates,N_cytosines,seqErr_true_rep,bsEff_true_rep,bsbEff_true_rep)


            #Save input for Stan when all cytosines are run together
            #Notice that the gamma distribution has different parameterization in Python and in Stan

            ind_R_p1=Z_R+1
            ind_C_p1=Z_C+1

            luxusv3_data={'n_cytosines': N_cytosines,'n_replicates': N_replicates*2,'n_predictors':N_predictors,
                'bsEff': bsEff_rep, 'bsBEff': bsbEff_rep, 'seqErr': seqErr_rep,
                'bsC': counts.astype(int), 'bsTot': reads_rep.astype(int),
                'X': numpy.tile(D,(N_cytosines,1)), 'Z_R': ind_R_p1.astype(int), 'Z_C': ind_C_p1.astype(int), 'alpha': a_E, 'beta': b_E,
                'alphaR': a_R, 'betaR': b_R,
                'alphal': 38, 'betal': 1,
                'sigmaB2': options.sigmaB2_stan, 'alphaC': a_C, 'betaC': b_C,
                'coordinates': coordinates[::N_predictors]}


            if options.save_LuxUS==1:
                output10=open("%s/LuxUS_simulated_all_%s_set%s_diff%s.txt"%(outputFolder,fileID,i,diff),'ab+')
                pickle.dump(luxusv3_data,output10)
                output10.close()


            #Form proportion table and design matrix
			

            #The proportion	table header preparation
            label1='control'
            label2='case'
            labels=list()
			
            indl1=1
            indl2=1


            for k1 in range(0,2*N_replicates):
                if D[k1,1]==1:
                    labels.append("%s%d"%(label2,indl1))
                    indl1+=1
                else:
                    labels.append("%s%d"%(label1,indl2))
                    indl2+=1

            #Add a newline in the end of the list 	
            labels.append('\n')
			

            coordinates_pt=coordinates[::N_predictors]
            if options.save_proportion_table==1:
                with open("%s/proportion_table_%s_set%s_diff%s.txt"%(outputFolder,fileID,i,diff),'a+') as pt:

                    #Write header line
                    pt.write('\t'.join(map(str,labels)))
				  
                    #Write data for each cytosine line by line
                    for k in range(0,N_cytosines): 

                        pt_row=list()
                        pt_row.append("chr1:%d:+:CpG"%(coordinates_pt[k]))
  
                        for k1 in range(0,N_replicates*2):
                            pt_row.append(N_reads)
                            pt_row.append(counts[k*N_replicates*2+k1].astype(int))

                        pt_row.append('\n')
                        pt.write('\t'.join(str(v) for v in pt_row))			

                with open("%s/design_matrix_%s_set%s_diff%s.txt"%(outputFolder,fileID,i,diff),'a+') as dm:
  
                    #Write header
                    dm.write("%s\t%s\n"%("Intercept","IsCase"))

                    #Write experimental design for each replicate separately
                    for k in range(0,N_replicates*2):
                        dm.write("%s\t%d\t%d\n"%(str(labels[k]),1,D[k,1].astype(int)))


            #Save input for LuxUS analysis where each cytosine is analysed separately

            if options.save_LuxUS_sep==1:
                for k in range(0,N_cytosines):
	
                    counts_sep=counts[range(k*N_replicates*2,(k+1)*N_replicates*2)]

                    luxsc_data_sep = {'n_replicates': N_replicates*2,'n_predictors':N_predictors,
                        'bsEff': bsEff_rep_sep, 'bsBEff': bsbEff_rep_sep, 'seqErr': seqErr_rep_sep,
                        'bsC': counts_sep.astype(int), 'bsTot': reads_rep_sep.astype(int),
                        'X': D, 'alpha': a_E, 'beta': b_E,
                        'sigmaB2': options.sigmaB2_stan, 'alphaR': a_R, 'betaR': b_R}

                    #Saving the generated datasets
                    output2=open("%s/LuxUS_simulated_sep_%s_cyt%s_set%s_diff%s.txt"%(outputFolder,fileID,k,i,diff),'ab+')
                    pickle.dump(luxsc_data_sep,output2)
                    output2.close()
