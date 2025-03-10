# Reducing GO terms visualization for heat response genes
library(AnnotationHub)
library(readxl)
library(GOSemSim)
library(rrvgo)
library(topGO)
library(AnnotationForge)

# Get A. millepora database: Technical reference: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html
hub <- AnnotationHub()
q <- query(hub, "Acropora millepora")
id <- q$ah_id[length(q)]
Amillepora <- hub[[id]]
# What key types are in the A. millepora database?
keytypes(Amillepora)

# "ACCNUM": GenBank accession numbers
# "ALIAS": Commonly used gene symbols
# "ENTREZID": Entrez gene Identifiers
# "GENENAME": The full gene name
# "GID": NA
# "PMID": Pubmed Identifiers
# "REFSEQ": Refseq Identifiers
# "SYMBOL": The official gene symbol

# Add the GO:ID to the database
geneID2GO <- topGO::readMappings(file="~/Lab Notebook/Chapter2/Data_analysis/Amil_GOmap.txt",sep="\t",IDsep=",")

# Look to see the columns 
head(select(Amillepora,
            keys = keys(Amillepora, keytype='GID'),
            keytype = 'GID',
            columns = columns(Amillepora)), n = 20)

# Load file that can translate 'Amillepora' gene names to a key in the ObjDB and to attach GO terms to the ObjDB, Amillepora


# GO Term similarity 
AmilleporaGO <- godata(Amillepora, ont="BP", keytype = "REFSEQ")
heat.response.GO.BP <- read.delim("~/Lab Notebook/Chapter2/Data_analysis/topGO_output/amillepora.hr.genes.GO_BP.csv", sep =',')
simMatrix <- calculateSimMatrix(heat.response.GO.BP$GO.ID,
                                orgdb = orgdb,
                                keytype = "GID",
                                semdata = GOSemSim::godata(orgdb, ont = "BP", keytype = "GID"),
                                ont="BP",
                                method="Rel")
