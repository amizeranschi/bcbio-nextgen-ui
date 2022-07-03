import json
import sys
import os.path

workflow_name = sys.argv[1]
path_to_images = sys.argv[2]
print(workflow_name)
print(path_to_images)

plot_list_GO = []
plot_list_KEGG = []
plot_list_disease = []
plot_list_mesh = []
go_analysis = False
kegg_analysis = False
disease_analysis = False
mesh_analysis = False
plot_type_list = []

if workflow_name == "variant_calling":
    # gene ontology
    if os.path.isfile(path_to_images + '/variant_calling_ego_high.jpg'):
        desc = "Enrichment Analysis - representation of gene ontology on HIGH impact genes"
        plot_list_GO += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_ego_high.jpg'}]
        go_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_ego_moderate.jpg'):
        desc = "Enrichment Analysis - representation of gene ontology on MODERATE impact genes"
        plot_list_GO += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_ego_moderate.jpg'}]
        go_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_ego_moderate_high.jpg'):
        desc = "Enrichment Analysis - representation of gene ontology on HIGH and MODERATE impact genes"
        plot_list_GO += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_ego_moderate_high.jpg'}]
        go_analysis =  True

    # kegg pathway
    if os.path.isfile(path_to_images + '/variant_calling_ego_high.jpg'):
        desc = "KEGG pathway over-representation analysis on the HIGH impact genes"
        plot_list_KEGG += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_kegg_high.jpg'}]
        kegg_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_kegg_moderate.jpg'):
        desc = "KEGG pathway over-representation analysis on the MODERATE impact genes"
        plot_list_KEGG += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_kegg_moderate.jpg'}]
        kegg_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_kegg_high_moderate.jpg'):
        desc = "KEGG pathway over-representation analysis on the HIGH and MODERATE impact genes"
        plot_list_KEGG += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_kegg_high_moderate.jpg'}]
        kegg_analysis =  True

    #  disease ontology
    if os.path.isfile(path_to_images + '/variant_calling_x_high.jpg'):
        desc = "Over-representation analysis for disease ontology in HIGH impact genes"
        plot_list_disease += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_x_high.jpg'}]
        disease_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_x_moderate.jpg'):
        desc = "Over-representation analysis for disease ontology in MODERATE impact genes"
        plot_list_disease += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_x_moderate.jpg'}]
        disease_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_x_moderate_high.jpg'):
        desc = "Over-representation analysis for disease ontology in HIGH and MODERATE impact genes"
        plot_list_disease += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_x_moderate_high.jpg'}]
        disease_analysis =  True

    #  mesh enrichment analysis
    if os.path.isfile(path_to_images + '/variant_calling_mesh_high.jpg'):
        desc = "MeSH enrichment analysis in HIGH impact genes"
        plot_list_mesh += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_mesh_high.jpg'}]
        mesh_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_mesh_moderate.jpg'):
        desc = "MeSH enrichment analysis in MODERATE impact genes"
        plot_list_mesh += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_mesh_moderate.jpg'}]
        mesh_analysis =  True
    if os.path.isfile(path_to_images + '/variant_calling_mesh_high_moderate.jpg'):
        desc = "MeSH enrichment analysis in HIGH and MODERATE impact genes"
        plot_list_mesh += [{'description': desc, 'imagePath': path_to_images + '/variant_calling_mesh_high_moderate.jpg'}]
        mesh_analysis =  True



    if go_analysis:
        plot_type_list += [{'title': "Gene Ontology Enrichment Analysis", 'components': plot_list_GO}]
    if kegg_analysis:
        plot_type_list += [{'title': "KEGG pathway over-representation analysis", 'components': plot_list_KEGG}]
    if disease_analysis:
        plot_type_list += [{'title': "Over-representation analysis for disease ontology", 'components': plot_list_disease}]
    if mesh_analysis:
        plot_type_list += [{'title': "MeSH enrichment analysis", 'components': plot_list_mesh}]

if workflow_name == "bulk_rna_seq":
    ma_analysis = False
    sample_analysis = False
    pca_analysis = False
    top_analysis = False
    plot_list_ma = []
    plot_list_heatmap = []
    plot_list_pca = []
    plot_list_top = []

    # MA plot
    if os.path.isfile(path_to_images + '/bulk_rna_seq_MA.png'):
        desc = "MA plot"
        plot_list_ma += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_MA.png'}]
        ma_analysis =  True
    # sample distances
    if os.path.isfile(path_to_images + '/bulk_rna_seq_heatmap_distances.png'):
        desc = "Heatmap of sample distances"
        plot_list_heatmap += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_heatmap_distances.png'}]
        sample_analysis =  True

    # PCA 
    if os.path.isfile(path_to_images + '/bulk_rna_seq_PCA-plot.png'):
        desc = "PCA plot"
        plot_list_pca += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_PCA-plot.png'}]
        pca_analysis =  True

    # top genes
    if os.path.isfile(path_to_images + '/bulk_rna_seq_top_genes.png'):
        desc = "Top significantly differentially expressed pheatmap plot"
        plot_list_top += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_top_genes.png'}]
        top_analysis =  True

    # gene ontology
    if os.path.isfile(path_to_images + '/bulk_rna_seq_ego_RT.png'):
        desc = "Enrichment Analysis - representation of gene ontology on ROOT genes"
        plot_list_GO += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_ego_RT.png'}]
        go_analysis =  True
    if os.path.isfile(path_to_images + '/bulk_rna_seq_ego_AE.png'):
        desc = "Enrichment Analysis - representation of gene ontology on AERIAL genes"
        plot_list_GO += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_ego_AE.png'}]
        go_analysis =  True

    # kegg pathway
    if os.path.isfile(path_to_images + '/bulk_rna_seq_kegg_RT.png'):
        desc = "KEGG pathway over-representation analysis on ROOT genes"
        plot_list_KEGG += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_kegg_RT.png'}]
        kegg_analysis =  True
    if os.path.isfile(path_to_images + '/bulk_rna_seq_kegg_AE.png'):
        desc = "KEGG pathway over-representation analysis on AERIAL genes"
        plot_list_KEGG += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_kegg_AE.png'}]
        kegg_analysis =  True

    #  disease ontology
    if os.path.isfile(path_to_images + '/bulk_rna_seq_x_RT.png'):
        desc = "Over-representation analysis for disease ontology in ROOT genes"
        plot_list_disease += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_x_RT.png'}]
        disease_analysis =  True
    if os.path.isfile(path_to_images + '/bulk_rna_seq_x_AE.png'):
        desc = "Over-representation analysis for disease ontology AERIAL genes"
        plot_list_disease += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_x_AE.png'}]
        disease_analysis =  True

    #  mesh enrichment analysis
    if os.path.isfile(path_to_images + '/bulk_rna_seq_mesh_RT.png'):
        desc = "MeSH enrichment analysis in ROOT genes"
        plot_list_mesh += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_mesh_RT.png'}]
        mesh_analysis =  True
    if os.path.isfile(path_to_images + '/bulk_rna_seq_mesh_AE.png'):
        desc = "MeSH enrichment analysis in AERIAL genes"
        plot_list_mesh += [{'description': desc, 'imagePath': path_to_images + '/bulk_rna_seq_mesh_AE.png'}]
        mesh_analysis =  True

    if ma_analysis:
        plot_type_list += [{'title': "MA plot", 'components': plot_list_ma}]
    if sample_analysis:
        plot_type_list += [{'title': "Heatmap of sample distances", 'components': plot_list_heatmap}]
    if pca_analysis:
        plot_type_list += [{'title': "PCA plot", 'components': plot_list_pca}]
    if top_analysis:
        plot_type_list += [{'title': "Top significantly differentially expressed pheatmap plot", 'components': plot_list_top}]
    if go_analysis:
        plot_type_list += [{'title': "Gene Ontology Enrichment Analysis", 'components': plot_list_GO}]
    if kegg_analysis:
        plot_type_list += [{'title': "KEGG pathway over-representation analysis", 'components': plot_list_KEGG}]
    if disease_analysis:
        plot_type_list += [{'title': "Over-representation analysis for disease ontology", 'components': plot_list_disease}]
    if mesh_analysis:
        plot_type_list += [{'title': "MeSH enrichment analysis", 'components': plot_list_mesh}]
if workflow_name == "atac_seq":
    plot_list = [{'description': "analysis_name", 'imagePath': path_to_images}]
    plot_type_list = [{'title': "analysis_name", 'components': plot_list}]


fh = open("report_context.json", "a+")
fh.write(json.dumps({"workflow": workflow_name, "content": plot_type_list}))
fh.close()