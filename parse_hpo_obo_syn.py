import pandas as pd
from collections import defaultdict
import csv
import sys

HPO_FILE = 'Archive/hp.obo.txt'
HPO_WITHOUT_SUPERCLASSES_FILE = 'hpo_without_superclasses.txt'
HPO_SUPERCLASSES_FILE = 'hpo_superclasses.txt'
HPO_NAMES_FILE = 'hpo_names.txt'
ALL_HPO_TERMS_AND_SYNONYMS_FILE = 'all_hpo_terms_and_synonyms.txt'
HPO_SUBCLASSES_FILE = 'hpo_subclasses.txt'

OLD_GENE_TO_PHENOTYPE_FILE = 'Archive/OMIM_FREQUENT_FEATURES_genes_to_phenotype.txt'
OLD_PHENOTYPE_TO_GENE_FILE = 'Archive/OMIM_FREQUENT_FEATURES_phenotype_to_genes.txt'
NEW_GENE_TO_PHENOTYPE_FILE = 'Expanded_OMIM_FREQUENT_FEATURES_genes_to_phenotype.txt'
NEW_PHENOTYPE_TO_GENE_FILE = 'Expanded_OMIM_FREQUENT_FEATURES_phenotype_to_genes.txt'
NEW_PHENOTYPE_TO_GENE_WITHOUT_SYNONYM_FILE = 'Expanded_OMIM_FREQUENT_FEATURES_phenotype_to_genes_without_synonym.txt'

'''
OLD_GENE_TO_PHENOTYPE_FILE = 'Archive/OMIM_ALL_FREQUENCIES_genes_to_phenotype.txt'
OLD_PHENOTYPE_TO_GENE_FILE = 'Archive/OMIM_ALL_FREQUENCIES_phenotype_to_genes.txt'
NEW_GENE_TO_PHENOTYPE_FILE = 'Expanded_OMIM_ALL_FREQUENCIES_genes_to_phenotype.txt'
NEW_PHENOTYPE_TO_GENE_FILE = 'Expanded_OMIM_ALL_FREQUENCIES_phenotype_to_genes.txt'
NEW_PHENOTYPE_TO_GENE_WITHOUT_SYNONYM_FILE = 'Expanded_OMIM_ALL_FREQUENCIES_phenotype_to_genes_without_synonym.txt'
'''

## Parse the original hpo file
res = defaultdict(dict)
with open(HPO_FILE, 'rb') as infile:
    START_OF_NEW_TERM = False  
    for line in infile.readlines():
        line = line.rstrip('\n\r')
        if line == '[Term]':
            START_OF_NEW_TERM = True
            continue
        if not line:
            START_OF_NEW_TERM = False
            continue
        if START_OF_NEW_TERM:
            parts = line.split(': ')
            key = parts[0]
            value = ': '.join(parts[1:]) 
            if key == 'id':
                print value 
                current_term = value
            if key == 'synonym':
                value = value.split('" ')[0][1:] 
            if key == 'is_a':
                superclass = value.split(' ! ')[0]
                if superclass in res['subclasses']:
                    res['subclasses'][superclass].append(current_term)
                else:
                    res['subclasses'][superclass] = [current_term]

            # keys may include ['comment', 'subset', 'xref', 'synonym', 'name', 'is_anonymous', 
            #                    'is_obsolete', 'alt_id', 'created_by', 'creation_date', 'is_a', 
            #                     'property_value', 'replaced_by', 'id', 'def', 'consider']

            if key in (['id', 'name', 'def', 'comment', 'subset', 'is_anonymous', 
                        'is_obsolete', 'created_by', 'creation_date', 'replaced_by', 'consider']):
                res[key][current_term] = value
            else:
                if current_term in res[key]:
                    res[key][current_term].append(value)
                else:
                    res[key][current_term] = [value]

keys = (['name', 'alt_id', 'def', 'comment', 'synonym', 'is_a', 'subclasses', 'xref', 'subset'])
dfs = []
for key in keys:
    df = pd.DataFrame(res[key].items(), columns = ['HPOID', key])
    dfs.append(df)
                    
df_res = pd.DataFrame()
for i in xrange(len(dfs)):
    temp_df = dfs[i]
    print temp_df.head()
    if df_res.empty:
        df_res = temp_df
    else:
        df_res = df_res.merge(temp_df, on = 'HPOID', how = 'outer')

hpo_without_superclasses =  df_res.loc[df_res.is_a.isnull(), ['HPOID', 'name']]
hpo_without_superclasses.to_csv(HPO_WITHOUT_SUPERCLASSES_FILE, sep = '\t', index = False)

superclasses = res['is_a']
for key in superclasses:
    value = superclasses[key]
    value = [_.split(' ! ')[0] for _ in value]
    superclasses[key] = value 
df_superclasses = pd.DataFrame(superclasses.items(), columns = ['hpoid', 'superclasses'])
df_superclasses.to_csv(HPO_SUPERCLASSES_FILE, sep = '\t', index = False)

hponames = res['name']
df_hponames = pd.DataFrame(hponames.items(), columns = ['hpoid', 'name'])
df_hponames.to_csv(HPO_NAMES_FILE, sep = '\t', index = False)


'''
Sample lines of df_res:

            HPOID          id                                               name alt_id                                                def comment  \
12296  HP:0005599  HP:0005599                           Hypopigmentation of hair    NaN                                                NaN     NaN   
12297  HP:0005598  HP:0005598  Facial telangiectasia in butterfly midface dis...    NaN  "Telangiectases (small dilated blood vessels) ...     NaN   
12298  HP:0030792  HP:0030792                                       Jaw neoplasm    NaN  "A tumor originating in the jaw (mandible or m...     NaN   

                                 synonym                                               is_a                            subclasses  \
12296            [Hair hypopigmentation]    [HP:0009887 ! Abnormality of hair pigmentation]  [HP:0001022, HP:0011358, HP:0011365]   
12297  [Butterfly facial telangiectasia]               [HP:0007380 ! Facial telangiectasia]                                   NaN   
12298                                NaN  [HP:0012289 ! Facial neoplasm, HP:0030791 ! Ab...                                   NaN   

                                                    xref        subset  
12296                                    [UMLS:C3278401]  hposlim_core  
12297                                    [UMLS:C4021632]           NaN  
12298  [MSH:D007573, SNOMEDCT_US:126634001, UMLS:C002...           NaN 
'''

## Export all hpo terms and their synonyms to a text file
tmp = df_res[['HPOID', 'name', 'synonym']]
hpo_res = []
for index, row in tmp.iterrows():
    hpoid = row['HPOID']
    name = row['name']
    synonym = row['synonym']
    hpo_res.append([hpoid, name])
    if type(synonym) == list:
        for syn in synonym:
            hpo_res.append([hpoid + '-synonym', syn])
df_hpo_res = pd.DataFrame(hpo_res, columns = ['hpoid', 'hponame']) 
df_hpo_res.sort(['hponame', 'hpoid'], ascending = [1, 1], inplace = True)
df_hpo_res.drop_duplicates(subset = ['hponame'], take_last = False, inplace = True)
df_hpo_res.to_csv(ALL_HPO_TERMS_AND_SYNONYMS_FILE, sep = '\t', index = False)

## Expand OMIM_ALL_FREQUENCIES_genes_to_phenotype.txt by associating the subclasses of each phenotype to the gene
def getAllChildren(node, tree):
    if node not in tree:
        return
    for child in tree[node]:
        yield child
        for grandchild in getAllChildren(child, tree):
            yield grandchild

tree = res['subclasses']
df_tree = pd.DataFrame(tree.items(), columns = ['hpoid', 'subclasses'])
df_tree.to_csv(HPO_SUBCLASSES_FILE, sep = '\t', index = False)
hpoid_name = res['name']
hpoid_synonym = res['synonym']

subclasses = dict()
geneid_to_symbol = dict()
expanded_hpoid2gene = []
with open(OLD_GENE_TO_PHENOTYPE_FILE, 'rb') as infile:
    f = csv.reader(infile, delimiter = '\t')
    firstline = True
    for line in f:
        if firstline:
            firstline = False
            continue
        geneid = line[0]
        genesymbol = line[1]
        hponame = line[2]
        hpoid = line[3]
        geneid_to_symbol[geneid] = genesymbol
        expanded_hpoid2gene.append([geneid, genesymbol, hponame, hpoid])
        # Add synonyms for this hpo term
        if hpoid in hpoid_synonym:
            for synonym in hpoid_synonym[hpoid]:
                expanded_hpoid2gene.append([geneid, genesymbol, synonym, hpoid + '-synonym'])

        # Get subclasses for this hpo term 
        try:
            subclass = list(getAllChildren(hpoid, tree))
            if (geneid, hpoid) in subclasses:
                subclasses[(geneid, hpoid)] += subclass
            else:
                subclasses[(geneid, hpoid)] = subclass
        except TypeError as e:
            print e
            pass

    for (geneid, hpoid) in subclasses:
        subclass = subclasses[(geneid, hpoid)]
        subclass = set(subclass)
        genesymbol = geneid_to_symbol[geneid]
        #if genesymbol == 'PTPN11':
        #    print hpoid, subclass
        for id in subclass:
            hponame = hpoid_name[id] 
            expanded_hpoid2gene.append([geneid, genesymbol, hponame, id])
            if id in hpoid_synonym:
                for synonym in hpoid_synonym[id]:              
                    expanded_hpoid2gene.append([geneid, genesymbol, synonym, id + '-synonym'])
 
df_expanded_hpoid2gene = pd.DataFrame(expanded_hpoid2gene, columns = ['geneid', 'genesymbol', 'hponame', 'hpoid']) 
df_expanded_hpoid2gene.sort(['genesymbol', 'hponame', 'hpoid'], ascending = [1, 1, 1], inplace = True)
df_expanded_hpoid2gene.drop_duplicates(subset = ['geneid', 'genesymbol', 'hponame'], take_last = False, inplace = True)
df_expanded_hpoid2gene.to_csv(NEW_GENE_TO_PHENOTYPE_FILE, sep = '\t', index = False)

# Expand the phenotype to gene file
df_old = pd.read_csv(OLD_PHENOTYPE_TO_GENE_FILE, sep = '\t', names = ['hpoid', 'hponame', 'geneid', 'genesymbol'], dtype = str)
df_old = df_old[1:]
df_old = df_old[['geneid', 'genesymbol', 'hponame', 'hpoid']]
df = pd.concat([df_old, df_expanded_hpoid2gene])
df.sort(['hponame', 'genesymbol', 'hpoid'], ascending = [1, 1, 1], inplace = True)
df.drop_duplicates(subset = ['geneid', 'genesymbol', 'hponame'], take_last = False, inplace = True)
df.to_csv(NEW_PHENOTYPE_TO_GENE_FILE, index = False, sep = '\t')

df_without_synonym = df[~df.hpoid.str.contains('-synonym')] 
df_without_synonym.to_csv(NEW_PHENOTYPE_TO_GENE_WITHOUT_SYNONYM_FILE, index = False, sep = '\t')
