

COUNTRY=CHAD STEPS="download filter samplesheet" bash submit_assemblies.sh
COUNTRY=Cameroon SPECIMEN_COL=specimen_type_x STEPS="download filter samplesheet" bash submit_assemblies.sh
COUNTRY=Morocco SPECIMEN_COL=specimen_type STEPS="download filter samplesheet " bash submit_assemblies.sh
COUNTRY=Zimbabwe COUNTRY_COL="Admin1" SPECIMEN_COL=specimen_type STEPS="download filter samplesheet" bash submit_assemblies.sh
COUNTRY=Ethiopia STEPS="download filter samplesheet" bash submit_assemblies.sh



chmod -R o+rX,go-w $DATA_DIR

COUNTRY=Morocco STEPS="run" bash submit_assemblies.sh
sbatch --export=ALL,COUNTRY=Ethiopia,STEPS=run submit_assemblies.sh 