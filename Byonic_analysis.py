###PYTHON#############
#Purpose : Build a DB of mutations found in the MS-tryspin/chemotrypsin digests
#Author : Aditya Ambati 
#date : May 26 2017
#update : version 1
#####import libraries
from itertools import chain
import pandas as pd
import numpy as np
%matplotlib
import re
from Bio import Entrez
from Bio import SeqIO


######################make a list of files 
arp_tryp = ['AFLPA319BB_ARP_TRYPSIN.csv','AFLPA328AA_ARP_TRYPSIN.csv','AFLPA359AA_ARP_TRYPSIN.csv','AFLPA373BA_ARP_TRYPSIN.csv','SF1B045CL_ARP_TRYPSIN.csv']
arp_chymo = ['AFLPA319BB_ARP_CHYMO.csv','AFLPA328AA_ARP_CHYMO.csv','AFLPA359AA_ARP_CHYMO.csv','AFLPA373BA_ARP_CHYMO.csv','SF1B0454CL_ARP_CHYMO.csv']
pan = ['AFLSA096AA_PAN_TRYPSIN.csv','AFLSA097AA_PAN_TRYPSIN.csv','AFLSA167AB_PAN_TRYPSIN.csv','AFLSA208A_PAN_TRYPSIN.csv', 'AFLSFDA280_PAN_TRYPSIN.csv', 'AFLSA174AA_PAN_TRYPSIN.csv']
pan_chymo = ['AFLSA208A_PAN_CHYMO.csv', 'AFLSFDA280_PAN_CHYMO.csv', 'AFLSA174AA_PAN_CHYMO.csv']
######### make a list of references influenza proteins used in the MS analysis
###############HA 				NP 			NA 				M1			PB2
references =['ADE29095.1', 'ADE29096.1', 'ADE29097.1', 'ADE29098.1', 'ADE29092.1'] ## to add in more id of references 
### define helper functions
def reference_retreive(references): ## this function will take alist and output a list with Dataframes for each individual peptide in a protein
	reference_proteins =[]
	for i in references:
		Entrez.email = "ambati@stanford.edu"
		handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=i)
		seq_record = SeqIO.read(handle, "fasta")
		tempref=pd.DataFrame(np.arange(1, len(seq_record.seq.tostring())+1))
		tempref.loc[:,('refseq')]=list(seq_record.seq.tostring())
		tempref.columns = ['position', 'refseq']
		tempref.loc[:, ('strain_protein')] = seq_record.description
		reference_proteins.append(tempref)
	return(reference_proteins)
################ this function will take a list with filenames and return a MS output spectra file with trimmed columns and calculated peptide length

def peponly(pepseq):
	return(re.sub(pattern='[^A-Z]*',repl='', string=pepseq))
def spectra_processing(spectra):
	dfs =[]
	for i in range(0, len(spectra)):
		
		temp = pd.read_csv(spectra[i], header = 0, sep = ";")
		temp1 = temp[['Protein\nRank','Peptide\n< ProteinMetrics Confidential >','Modification Type(s)', 'Starting\nposition', 'Protein Name', '|Log Prob|']]

	## rename the columns
		temp1.columns = ['protein_rank', 'peptide', 'modtype', 'start_pos', 'protein_name', 'log_prob']
		temp2=temp1.loc[~temp1['protein_name'].str.contains('Reverse')]
		temp2.reset_index(drop=True, inplace=True)

		temp2.loc[:, ('filename')] = spectra[i]
		temp2.loc[:, ('pep_wt_tags')]=[re.findall('\..*\.', i) for i in temp2.peptide] ## take only characters between the first and the last dot character
		temp2.loc[:, ('tag_rem_pep')] = [peponly(i) for i in temp2.pep_wt_tags] ### assign to a new column the raw peptide without tags
		temp2.loc[:, ('tag_rm_len')] = [len(i) for i in temp2.tag_rem_pep] ## get the raw peptide length
		temp2.loc[:, ('templen')] = temp2.tag_rm_len + temp2.start_pos ### add the startingpos + untagged peptide length
		#temp2.loc[:, 'org_length'] = [len(i) for i in temp2.peptide] 
		dfs.append(temp2)
	return(dfs)
################################################### this function will update the reference files 

def db_reference(df, reference, i): #### this is an old function, use the one below called mutcheck!
	i = int(i)
	index_df = np.arange(df.start_pos[i]-1, df.templen[i]-1) ## get index from parent for loop
	reference.loc[index_df, 'total_occurences'] = reference.loc[index_df, 'total_occurences']+1 # if index from spectra is present in reference index add 1 each time to all index positions
	ref_index = reference.ix[index_df, ('refseq')] ## get the reference index peptides
	df_pep = pd.Series(list(df.tag_rem_pep[i]), index=ref_index.index) ### get the spectra index peptides
	ref_mutindex = ref_index.index[ref_index != list(df.tag_rem_pep[i])] ## if df spectra index peptides are not equal to the reference index peptides, get that ref_mut index
	reference.loc[ref_mutindex, ('mutated_occurence')] = reference.loc[ref_mutindex, ('mutated_occurence')] +1 ## if there appears to be a mutation add 1 each time it finds one
	reference.loc[ref_mutindex, ('mutated_AA')] = df_pep[df_pep != ref_index] ### finally retreive the mutated pos from spectra index peptide and add it into mutated pos in the reference file

####################
def mutcheck(spectraDF, reference, stringquery):
	tempdf = spectraDF[spectraDF['protein_name'].str.contains(stringquery)]
	tempdf.reset_index(drop=False, inplace =True)
	reference.loc[:, 'mutated_AA'] = 'Nan'
	reference.loc[:, 'total_occurences'] =0
	reference.loc[:, 'mutated_occurence']=0
	reference.loc[:, 'occu_log_prob'] =0
	reference.loc[:, 'mut_log_prob']=0
	reference.loc[:, 'fileID'] = tempdf.filename.unique()
	mutcheck =[]
	for i in xrange(0, len(tempdf.start_pos)):
		HAindex = np.arange(tempdf.start_pos[i]-1, tempdf.templen[i]-1)
		reference.loc[HAindex, 'total_occurences'] = reference.loc[HAindex, 'total_occurences'] +1
		reference.loc[HAindex, 'occu_log_prob'] = reference.loc[HAindex, 'occu_log_prob'] + tempdf.log_prob[i]
		q=reference.ix[HAindex, ('refseq')] 
		MSpep = pd.Series(list(tempdf.tag_rem_pep[i]), index=q.index)
		mutindex=q.index[q != list(tempdf.tag_rem_pep[i])]
		if len(mutindex) > 1:
			mutcheck.append(i)
		reference.loc[mutindex, ('mutated_occurence')] = reference.loc[mutindex, ('mutated_occurence')] +1
		reference.loc[mutindex, ('mut_log_prob')] = reference.loc[mutindex, ('mut_log_prob')] + tempdf.log_prob[i]
		reference.loc[mutindex, ('mutated_AA')] = MSpep[MSpep != q]

	return([tempdf.ix[mutcheck], reference])

##### process the spectra files 
pan_tryp_dfs = spectra_processing(pan)
pan_chymo_dfs= spectra_processing(pan_chymo)
arp_tryp_dfs = spectra_processing(arp_tryp)
arp_chymo_dfs = spectra_processing(arp_chymo)
######retreive these protein sequences 
ref_proteins = reference_retreive(references) ### feed a list of protein IDs to retreive them from database amd get them in a database format



###################################################################
HA_ref = ref_proteins[0][:]
NP_ref = ref_proteins[1][:]
NA_ref = ref_proteins[2][:]
MA_ref = ref_proteins[3][:]
PB2_ref = ref_proteins[4][:]

PB2check=mutcheck(pand_all, reference = PB2_ref[:], stringquery ='>PB2-|ACQ63273.1')
HAcheck=mutcheck(pand_all, reference = HA_ref[:], stringquery = '>HA-|ACQ55359.1|ADE29085.1')
NAcheck=mutcheck(pand_all, reference = NA_ref[:], stringquery = '>NA-|ADE29087.1')
MAcheck=mutcheck(pand_all, reference = MA_ref[:], stringquery = '>MA1-|ACP44184.1')
NPcheck=mutcheck(pand_all, reference = NP_ref[:], stringquery = '>NP-|ACP44183.1')

pan_tryp_DB =[]
pan_checkids=[]
for df in pan_tryp_dfs:
	PB2check=mutcheck(df, reference = PB2_ref[:], stringquery ='>PB2-|ACQ63273.1')
	HAcheck=mutcheck(df, reference = HA_ref[:], stringquery = '>HA-|ACQ55359.1|ADE29085.1')
	NAcheck=mutcheck(df, reference = NA_ref[:], stringquery = '>NA-|ADE29087.1')
	MAcheck=mutcheck(df, reference = MA_ref[:], stringquery = '>MA1-|ACP44184.1')
	NPcheck=mutcheck(df, reference = NP_ref[:], stringquery = '>NP-|ACP44183.1')
	alldb=pd.concat([PB2check[1], HAcheck[1], NAcheck[1], MAcheck[1], NPcheck[1]], ignore_index=True)
	allcheck = pd.concat([PB2check[0], HAcheck[0], NAcheck[0], MAcheck[0], NPcheck[0]], ignore_index=True)

	pan_tryp_DB.append(alldb)
	pan_checkids.append(allcheck)

arp_tryp_DB=[]
arp_checkids=[]
for df in arp_tryp_dfs:
	PB2check=mutcheck(df, reference = PB2_ref[:], stringquery ='>PB2-|ACQ63273.1')
	HAcheck=mutcheck(df, reference = HA_ref[:], stringquery = '>HA-|ACQ55359.1|ADE29085.1')
	NAcheck=mutcheck(df, reference = NA_ref[:], stringquery = '>NA-|ADE29087.1')
	MAcheck=mutcheck(df, reference = MA_ref[:], stringquery = '>MA1-|ACP44184.1')
	NPcheck=mutcheck(df, reference = NP_ref[:], stringquery = '>NP-|ACP44183.1')
	alldb=pd.concat([PB2check[1], HAcheck[1], NAcheck[1], MAcheck[1], NPcheck[1]], ignore_index=True)
	allcheck = pd.concat([PB2check[0], HAcheck[0], NAcheck[0], MAcheck[0], NPcheck[0]], ignore_index=True)

	arp_tryp_DB.append(alldb)
	arp_checkids.append(allcheck)


arp_db=pd.concat(arp_tryp_DB, ignore_index=True)
pan_db=pd.concat(pan_tryp_DB, ignore_index=True)
arp_db.to_csv('/home/labcomp/Desktop/Vaccine_MS_May16/Arepanrix_trypsin_june16.csv')
pan_db.to_csv('/home/labcomp/Desktop/Vaccine_MS_May16/Pandemrix_trypsin_june16.csv')
arp_check=pd.concat(arp_checkids, ignore_index=True)
pan_check=pd.concat(pan_checkids, ignore_index=True)
arp_check
pan_check
pd.concat([arp_check, pan_check], ignore_index=True)
all_tryp_check=pd.concat([arp_check, pan_check], ignore_index=True)
all_tryp_check.to_csv('/home/labcomp/Desktop/Vaccine_MS_May16/check_back_excel raw_june16.csv')


############################################## post processing the DB files ########################################
import pandas as pd
import numpy as np

arp = pd.read_csv('Arepanrix_trypsin_june16.csv')
pan = pd.read_csv('Pandemrix_trypsin_june16.csv')
prot=pan.strain_protein.unique()
refprot=[i.split('[')[0] for i in prot] 
refid=[i.split(' ')[0] for i in refprot]
refprotname= ['PB2', 'HA', 'NA', 'MA1', 'NP'] 
### try to simplify the protein ids
['ADE29092.1 polymerase PB2 ',
 'ADE29095.1 hemagglutinin ',
 'ADE29097.1 neuraminidase ',
 'ADE29098.1 matrix protein 1 ',
 'ADE29096.1 nucleocapsid protein ']
#### change the long protein name to short ones for e.g NucleoProtein > NP etc.
arp.loc[:, 'prot'] = 'NaN'
pan.loc[:, 'prot'] = 'NaN'
for i, j in enumerate(refid): 
	arp.loc[arp.strain_protein.str.contains(j), ('prot')] = refprotname[i]
	pan.loc[pan.strain_protein.str.contains(j), ('prot')] = refprotname[i]


## drop uneceassary cols 
arp2=arp.drop(['Unnamed: 0', 'strain_protein'], axis=1)
pan2=pan.drop(['Unnamed: 0', 'strain_protein'], axis=1)
#### add keys column to these files by combining the proteinname and the position of that peptide in the ref protein
pan2.loc[:,('keys')]=[str(pan2.prot[i]) + '_' + str(pan2.position[i]) for i in xrange(0, len(pan2))]
arp2.loc[:,('keys')]=[str(arp2.prot[i]) + '_' + str(arp2.position[i]) for i in xrange(0, len(arp2))]

######## groupby pandas if a rows have the same common 'prot','position', 'refseq', 'mutated_AA' then mean the total_occurences', 'mutated_occurence', 'mut_log_prob'

pan2_grp2=pan2[['prot','position', 'refseq','mutated_AA','total_occurences', 'mutated_occurence', 'mut_log_prob']].groupby(['prot','position', 'refseq', 'mutated_AA']).mean()
pan2_reduced=pan2_grp2.reset_index()
#### further if certain rows satisfy conditions total_occurences > 0, mutated_occurence > = 3 and a mut_log prob of > 2, take these as valid rows 
pan2_final=pan2_reduced.loc[(pan2_reduced['total_occurences'] != 0) & (pan2_reduced['mutated_occurence'] >= 3) & (pan2_reduced['mut_log_prob'] > 2)]
pan2_sort=pan2_final.sort_values(['mutated_occurence'], ascending=False)
#pan2_sort.loc[:, ('percentage')]=(pan2_sort.mutated_occurence/pan2_sort.total_occurences)*100


## do the same with arpanrix
arp2_grp2=arp2[['prot','position', 'refseq','mutated_AA','total_occurences', 'mutated_occurence', 'mut_log_prob']].groupby(['prot','position', 'refseq', 'mutated_AA']).mean()
arp2_reduced=arp2_grp2.reset_index()
arp2_final=arp2_reduced.loc[(arp2_reduced['total_occurences'] != 0) & (arp2_reduced['mutated_occurence'] >= 3) & (arp2_reduced['mut_log_prob'] > 2)]
arp2_sort=arp2_final.sort_values(['mutated_occurence'], ascending=False)
#arp2_sort.loc[:, ('percentage')]=(arp2_sort.mutated_occurence/arp2_sort.total_occurences)*100

### reset index on these files but keep the column structure
pan2_sort.reset_index(drop=True, inplace=True)
arp2_sort.reset_index(drop=True, inplace=True)

#### get all key mutations from each of the vaccine for e.g HA 146 becomes HA_146
pan2_keys = [str(pan2_sort.prot[i]) + '_' + str(pan2_sort.position[i]) for i in xrange(0, len(pan2_sort))]
arp2_keys = [str(arp2_sort.prot[i]) + '_' + str(arp2_sort.position[i]) for i in xrange(0, len(arp2_sort))]
common_keys= set(pan2_keys + arp2_keys) #### from the selected rows in above sort files pan2_sort and arp2_sort, make a big list containing common unique keys
### make these keys from the orginal database file *_reduced.
pan2_reduced.loc[:, ('keys')] = [str(pan2_reduced.prot[i]) + '_' + str(pan2_reduced.position[i]) for i in xrange(0, len(pan2_reduced))]
arp2_reduced.loc[:, ('keys')] = [str(arp2_reduced.prot[i]) + '_' + str(arp2_reduced.position[i]) for i in xrange(0, len(arp2_reduced))]

###get all the positions from the main DB *_reduced present in the common keys
panf=pan2_reduced[pan2_reduced['keys'].isin(common_keys)]
arpf=arp2_reduced[arp2_reduced['keys'].isin(common_keys)]
#### here we find that certain mutations maybe different across same vaccine lots (for e.g within pandemrix at same position , batch1 might have a mutation at HA 146 as D whereas batch might have mutatation ar HA 146 as F)
### further, if this is the case this will showup in above at groupby as two rows because of the difference in mutated_AA. we are merely gathering all these empty keys and removing them
arpf2=arpf.index[(arpf.duplicated(['keys'], keep=False)) & (arpf.mutated_AA.str.contains('Nan'))]
panf2=panf.index[(panf.duplicated(['keys'], keep=False)) & (panf.mutated_AA.str.contains('Nan'))]
### drop these indexes from the mutated files arp_mut/pan_mut
arp_mut=arpf.drop(arpf2, axis=0)
pan_mut=panf.drop(panf2, axis=0)
#### merge the tow vaccines on keys and get suffixes for each of them.
all_mut=pd.merge(arp_mut, pan_mut, on='keys', how='outer', suffixes=['_ARP', '_PAN'])
#all_mut.to_csv('/Users/adityaambati/Desktop/All_mutations.csv')

#all_mut.loc[(all_mut.mutated_occurence_ARP >= 3) | (all_mut.mutated_occurence_PAN >=3)]
#all_mut[(all_mut.mutated_occurence_ARP >= 3) & (all_mut.mutated_occurence_PAN >=3)]

#all_mut.loc[(all_mut.mut_log_prob_ARP >= 2) & (all_mut.mut_log_prob_PAN >=2)]
##make a copy of the file
all_mut2 = all_mut
#### calculate the ratios mut/total *100 - gives %
all_mut2.loc[:, ('mutated_occurence_ARP')] = np.round(all_mut2.mutated_occurence_ARP)
all_mut2.loc[:, ('mutated_occurence_PAN')] = np.round(all_mut2.mutated_occurence_PAN)
all_mut2.loc[:, ('total_occurences_ARP')] = np.round(all_mut2.total_occurences_ARP)
all_mut2.loc[:, ('total_occurences_PAN')] = np.round(all_mut2.total_occurences_PAN)

all_mut2.loc[:, ('percent_ARP')]=np.round((all_mut2.mutated_occurence_ARP/all_mut2.total_occurences_ARP)*100, decimals=2)
all_mut2.loc[:, ('percent_PAN')]=np.round((all_mut2.mutated_occurence_PAN/all_mut2.total_occurences_PAN)*100, decimals=2)

### clean up cols
all_mut3=all_mut2.drop(['prot_PAN', 'position_PAN','mut_log_prob_ARP', 'mut_log_prob_PAN'], axis=1)
all_mut3.reset_index(drop=True, inplace=True)
all_mut3.loc[:, ('mut_pos_PAN')] = [all_mut3.refseq_PAN[i] + ' > ' + all_mut3.mutated_AA_PAN[i] for i in xrange(0, len(all_mut3))]
all_mut3.loc[:, ('mut_pos_ARP')] = [all_mut3.refseq_ARP[i] + ' > ' + all_mut3.mutated_AA_ARP[i] for i in xrange(0, len(all_mut3))]


#############Debugging !!!!!!!
### debugging and verifying mutations in the main DF raw files
pan_raw=pd.concat(pan_tryp_dfs, ignore_index=True)
arp_raw=pd.concat(arp_tryp_dfs, ignore_index=True)



arp_raw.loc[(arp_raw.protein_name.str.contains('>HA-|ACQ55359.1|ADE29085.1|hema')) & (arp_raw.start_pos == 137)]
pan_raw.loc[(pan_raw.protein_name.str.contains('>HA-|ACQ55359.1|ADE29085.1|hema')) & (pan_raw.start_pos == 137)]


### reach into the reference seq to get 20 aa around the mutations
HA_seq=''.join(ref_proteins[0].refseq.tolist())
NP_seq = ''.join(ref_proteins[1].refseq.tolist())
NA_seq = ''.join(ref_proteins[2].refseq.tolist())
MA_seq = ''.join(ref_proteins[3].refseq.tolist())
PB2_seq = ''.join(ref_proteins[4].refseq.tolist())

all_seq=pd.DataFrame([HA_seq, NP_seq, NA_seq, MA_seq, PB2_seq], columns=['seq'])
all_seq.loc[:, ('prots')] = ['HA', 'NP', 'NA', 'MA1', 'PB2']

all_mut3.loc[:, ('20AA_mut_ARP')] = 'NaN'
all_mut3.loc[:, ('20AA_mut_PAN')] = 'NaN'
all_mut3.loc[:, ('20AA_X179A')] = 'NaN'


for i in xrange(0, len(all_mut3)):
	pos = all_mut3.position_ARP[i]-1
	refprot = ''.join(all_seq.loc[all_seq.prots.str.contains(all_mut3.prot_ARP[i]), ('seq')].tolist())
	print refprot[pos-10:pos] + '['+all_mut3.mut_pos_ARP[i]+']' + refprot[pos+1:pos+10], refprot[pos-10:pos] + '['+all_mut3.mut_pos_PAN[i]+']' + refprot[pos+1:pos+10], refprot[pos-10:pos+10]
	all_mut3.loc[i, ('20AA_mut_ARP')] = refprot[pos-10:pos] + '['+all_mut3.mut_pos_ARP[i]+']' + refprot[pos+1:pos+10]
	all_mut3.loc[i, ('20AA_mut_PAN')] = refprot[pos-10:pos] + '['+all_mut3.mut_pos_PAN[i]+']' + refprot[pos+1:pos+10]
	all_mut3.loc[i, ('20AA_X179A')] = refprot[pos-10:pos+10]


################### go to orginal DB files and retreive raw counts
pan2_mut_raw_counts=pan2.loc[pan2['keys'].isin(common_keys)]
arp2_mut_raw_counts=arp2.loc[arp2['keys'].isin(common_keys)]


arp2_mut_raw_counts=arp2_mut_raw_counts.sort_values(['keys'])
pan2_mut_raw_counts=pan2_mut_raw_counts.sort_values(['keys'])

