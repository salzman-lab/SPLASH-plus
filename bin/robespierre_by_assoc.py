import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator

#### Take as input a compactor_classification script.
#### Take as input a compactor contingency table.
#### User needs only specify a donor and a tissue.

donor, tissue=sys.argv[1], sys.argv[2]
anchor_list_path='anchor.geq1tissue.geq6var'
anchor_list=pd.read_csv(anchor_list_path, engine='python',sep='\t')
path='/oak/stanford/groups/horence/Roozbeh/NOMAD/compactors/TSP_SS2_tissue/'
#classified='/oak/stanford/groups/horence/julias/V3/uni.some_compactors_aligned.tab'
#classified='/oak/stanford/groups/horence/julias/V3/anchor.geq1tissue.geq6var'
contingency1 = pd.read_csv(path + donor + '/' + tissue + '/sample_specificity.tsv', engine='python', sep='\t')
metadata_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/utility_files/meta_data/Tabula_Sapiens/TSP1_TSP15_metadata.csv"
metadata = pd.read_csv(metadata_file,sep=',',engine='python')
metadata['cell_id'] = [i.split('.')[0] for i in list(metadata['cell_id'])]
#metadata = metadata[metadata['cell_id'].str.contains('SS2')]
#classified='/oak/stanford/groups/horence/julias/V3/uni.some_compactors_aligned.tab'
classified=path+donor+'/'+tissue+'/'+'classified_compactors.tsv'
classified=pd.read_csv(classified, sep='\t',engine='python',quoting=3)

contingency_columns = [contingency1.columns[i].split('/')[-1].split('_R1_')[0] for i in range(len(contingency1.columns))]
contingency_columns = [contingency_columns[i].split('_R1.f')[0] for i in range(len(contingency_columns))]
# contingency_columns = contingency_columns[0:2] + contingency_columns[3:]
contingency1 = contingency1.iloc[:,:-1]
contingency1.columns = contingency_columns

contingency1 = contingency1.rename(columns={'compactor':'compactor_majority'})
mapper = classified[['compactor_valid','compactor_majority']]
contingency1 = contingency1.merge(mapper)#,how='left',on=['compactor_majority'])

contingency = contingency1
#classifiedd=classified[classified['n.var.person.anyspl']>4].fillna(0)
classifiedd=classified
#classifiedd = classifiedd[classifiedd['num_compactor_per_anchor']>=30]
#classifiedd = classifiedd[classifiedd['num_compactor_gene_anchor']==1]
classifiedd = classifiedd[classifiedd['run_length_A']<5]
classifiedd = classifiedd[classifiedd['run_length_C']<5]
classifiedd = classifiedd[classifiedd['run_length_G']<5]
classifiedd = classifiedd[classifiedd['run_length_T']<5]


#cand_genes = list(classified['cg'].dropna().unique())
#list_comp_lists = classifiedd.groupby(['anchor'])['compactor_valid'].apply(list)
trunc=donor+'_'+tissue+'_'
trunc=trunc[:10]
#list_comp_lists = list(anchor_list[anchor_list['myff']==trunc]['anchor'])
gene=sys.argv[3]
#list_comp_lists = list(classifiedd[classifiedd['compactor_gene'].str.contains(gene,case=False,na=False)]['anchor'].unique())
#into = classifiedd[classifiedd['compactor_gene'].str.contains(gene,case=False,na=False)]
anchors_take = anchor_list[anchor_list['compactor_gene'].str.contains(gene, case=False, na=False)]
anchors_take = anchors_take[anchors_take['myff'].str.contains(trunc, case=False, na=False)]
list_comp_lists = list(anchors_take['anchor'].unique())
#list_comp_lists = into.groupby(['compactor_gene'])['anchor'].apply(list)
#dna_colors = ['#ffffff','#8dd3c7','#ffffb3','#fb8072','#80b1d3',"lightsteelblue"] #["N","A","T","C","G"]
dna_colors = ['#ffffff','seagreen','#ffffb3','salmon','lightsteelblue']#,"lightsteelblue"] #["N","A","T","C","G"]

col_dict={
    -1:dna_colors[0],
    0:dna_colors[1],
    1:dna_colors[2],
    2:dna_colors[3],
    3:dna_colors[4]}#,
    #6:dna_colors[5]}

labels = np.array(["N","A","T","C","G"])#,"agree"])
len_lab = len(labels)
norm_bins = np.sort([*col_dict.keys()]) + 0.5
norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
diff = norm_bins[1:] - norm_bins[:-1]
tickz = norm_bins[:-1] + diff / 2
cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

def bpToInt(x):
    if x=='A':
        return 0
    if x=='T':
        return 1
    if x=='C':
        return 2
    if x=='G':
        return 3
    return 4

cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

for i in range(len(list_comp_lists)):
    matrix=contingency[contingency['anchor'] == list_comp_lists[i]]
    colmap=metadata[metadata['cell_id'].isin(matrix.columns)][['cell_id','cell_ontology_class']]
    ontologie=list(colmap['cell_ontology_class'].unique())
    if matrix.shape[0] >= 4:
        plottt = pd.DataFrame()
        sns.set(font_scale=1)
        for ont in ontologie:
            plottt[ont] = matrix[[i for i in matrix.columns if i in list(metadata[metadata['cell_ontology_class']==ont]['cell_id'])]].sum(axis=1)
        plottt['anchor'] = matrix['anchor']
        plottt['compactor_valid'] = matrix['compactor_valid']
        plottt = plottt.groupby(['anchor','compactor_valid']).sum()
        dummy = matrix.groupby(['anchor','compactor_valid']).sum()
        if dummy.shape[0] >= 4:
            """
            typefractions=plottt.div(plottt.sum(axis=0), axis=1)
            ax=sns.heatmap(typefractions,linewidths=0.5,annot=False,cmap='Blues')
            plt.title(donor + ' ' + tissue + ' ' + anchor_seq + ' ' + gene_name)
            plt.savefig(donor+'_'+tissue+'_'+gene_name+'_'+str(i)+'2.pdf',bbox_inches='tight')
            plt.close()
            """

            anchor_seq=str(list_comp_lists[i])
            #gene_name=str(list(classifiedd[classifiedd['anchor']==list_comp_lists[i][0]]['compactor_gene'])[0])
            #gene_name=str(list_comp_lists.index[i])
            gene_name=str(list(anchor_list[anchor_list['anchor']==list_comp_lists[i]]['compactor_gene'])[0])
            sns.set(font_scale=0.45)
            matrix=matrix.groupby(['anchor','compactor_valid']).sum()
            matrix=matrix.sort_values('compactor_valid',ascending=True)
            matrix=matrix.replace(0,np.nan).dropna(how='all', axis=1)
            try:
                matrix = matrix.sort_values(matrix.max().idxmax(), ascending=False)
            except ValueError:
                continue
            matrix = matrix.T
            matrix = matrix.sort_values(matrix.columns[0], ascending=False).T
            compactors = [i[1] for i in matrix.index]
            jar = classifiedd[classifiedd['compactor_valid'].isin(compactors)][['anchor','compactor_valid','compactor_gene']].set_index(['anchor','compactor_valid'])
            jar['STAR_annotated'] = [i!=0 for i in jar['compactor_gene'].fillna(0)]
            annot = jar['STAR_annotated']
            lut1 = dict(zip(annot.unique(), ['slateblue','firebrick']))
            row_colors1 = annot.map(lut1)
            ax=sns.clustermap(matrix.fillna(0),linewidths=0.5,annot=False,cmap='Blues', dendrogram_ratio=0.2,row_colors=row_colors1)
            ax.fig.set_size_inches(24,12)
            ax.ax_row_dendrogram.set_visible(False) #suppress row dendrogram
            ax.ax_col_dendrogram.set_visible(False) #suppress column dendrogram
            plt.title(donor + ' ' + tissue + ' ' + anchor_seq + ' ' + gene_name)
            plt.savefig(anchor_list_path+'.dir/single_cell/'+donor+'_'+tissue+'_'+gene_name+'_'+str(i)+'3.pdf',bbox_inches='tight')
            order = ax.dendrogram_row.reordered_ind
            ind = matrix.index[order]
            plt.close()

            compactor_list=[i[1] for i in ind]
            compactor_list = [i[27:] for i in compactor_list]
            n=max(len(x) for x in compactor_list)
            m=int(len(compactor_list))
            countsmat=-1*np.ones((m,n))
            for i,cons in enumerate(compactor_list):
                countsmat[i,:len(cons)] = [bpToInt(x) for x in cons]
            #for i in range(countsmat.shape[1]):
              #  if all(elem == countsmat[0,i] for elem in countsmat[:,i]):
               #     countsmat[:,i] = 6
            im = plt.imshow(countsmat,aspect='auto', norm=norm, cmap=cm,interpolation='nearest')
            plt.xlabel('sample')
            plt.ylabel('compactor sequence')
            plt.xlabel('sequence index')
            plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
            cb = plt.colorbar(im,orientation="horizontal",pad=0.2) ### just to get cb.ax
            cb = plt.colorbar(im, format=fmt, ticks=tickz, orientation="horizontal", cax=cb.ax)
            plt.gca().yaxis.tick_right()
            plt.yticks(ticks=range(len(compactor_list)),labels=compactor_list)
            plt.grid(b=None)
            plt.savefig(anchor_list_path+'.dir/single_cell/'+donor+'_'+tissue+'_'+gene_name+'_'+str(i)+'3_sequence.pdf',bbox_inches='tight')
            plt.close()

            """
            fractions=matrix.div(matrix.sum(axis=0), axis=1)
            ax=sns.heatmap(fractions,linewidths=0.5,annot=False,cmap='Blues')
            plt.title(donor + ' ' + tissue + ' ' + anchor_seq + ' ' + gene_name)
            plt.savefig(donor+'_'+tissue+'_'+gene_name+'_'+str(i)+'4.pdf',bbox_inches='tight')
            plt.close()
            """

            plottt = plottt.sort_values('compactor_valid', ascending=True)
            try:
                plottt = plottt.sort_values(plottt.max().idxmax(), ascending=False)
            except ValueError:
                continue
            plottt = plottt.T
            plottt = plottt.sort_values(plottt.columns[0], ascending=False).T
            plottt=plottt.replace(0,np.nan).dropna(how='all', axis=1)
            compactors = [i[1] for i in plottt.index]
            jar = classifiedd[classifiedd['compactor_valid'].isin(compactors)][['anchor','compactor_valid','compactor_gene']].set_index(['anchor','compactor_valid'])
            jar['STAR_annotated'] = [i!=0 for i in jar['compactor_gene'].fillna(0)]
            annot = jar['STAR_annotated']
            lut1 = dict(zip(annot.unique(), ['slateblue','firebrick']))
            row_colors1 = annot.map(lut1)
            ax=sns.clustermap(plottt.fillna(0),linewidths=0.5,annot=False,cmap='Blues', dendrogram_ratio=0.2,row_colors=row_colors1)
            ax.fig.set_size_inches(24,12)
            ax.ax_row_dendrogram.set_visible(False) #suppress row dendrogram
            ax.ax_col_dendrogram.set_visible(False) #suppress column dendrogram
            sns.set(rc={'figure.figsize':(11.27,11.27)})
            #x0, _y0, _w, _h = ax.cbar_pos
            #g.ax_cbar.set_position([1, 0.9, ax.ax_col_dendrogram.get_position().width, 0.01])
            plt.title(donor + ' ' + tissue + ' ' + anchor_seq + ' ' + gene_name)
            plt.savefig(anchor_list_path+'.dir/pseudobulked/'+donor+'_'+tissue+'_'+gene_name+'_'+str(i)+'1.pdf', bbox_inches='tight')
            order = ax.dendrogram_row.reordered_ind
            ind = plottt.index[order]
            plt.close()

            compactor_list=[i[1] for i in ind]
            compactor_list = [i[27:] for i in compactor_list]
            n=max(len(x) for x in compactor_list)
            m=int(len(compactor_list))
            countsmat=-1*np.ones((m,n))
            for i,cons in enumerate(compactor_list):
                countsmat[i,:len(cons)] = [bpToInt(x) for x in cons]
            #for i in range(countsmat.shape[1]):
             #   if all(elem == countsmat[0,i] for elem in countsmat[:,i]):
             #       countsmat[:,i] = 6
            im = plt.imshow(countsmat,aspect='auto', norm=norm, cmap=cm,interpolation='nearest')
            plt.xlabel('sample')
            plt.ylabel('compactor sequence')
            plt.xlabel('sequence index')
            plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
            plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
            cb = plt.colorbar(im,orientation="horizontal",pad=0.2) ### just to get cb.ax
            cb = plt.colorbar(im, format=fmt, ticks=tickz, orientation="horizontal", cax=cb.ax)
            plt.gca().yaxis.tick_right()
            plt.yticks(ticks=range(len(compactor_list)),labels=compactor_list)
            plt.grid(b=None)
            plt.savefig(anchor_list_path+'.dir/pseudobulked/'+donor+'_'+tissue+'_'+gene_name+'_'+str(i)+'1_sequence.pdf',bbox_inches='tight')
            plt.close()

            print('anchor '+str(i),flush=True)
