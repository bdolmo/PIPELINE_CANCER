# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/header.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Biomarkers.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/logos.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Other_variants.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Therapeutic_variants.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/rare_variants.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Summary_qc.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Sample_information.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Disclaimers.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/Disclaimers2.jrxml
# ./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/LungCancer_Report_v1.jrxml
sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/header.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/header.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Biomarkers_cat.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Biomarkers_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/logos.jrxml
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/logos.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Other_variants_cat.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Other_variants_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Therapeutic_variants_cat.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Therapeutic_variants_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/rare_variants_cat.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/rare_variants_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Summary_qc_cat.jrxml

./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Summary_qc_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g'  /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Sample_information_cat.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Sample_information_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Disclaimers_cat.jrxml
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Disclaimers_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Disclaimers2_cat.jrxml
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/Disclaimers2_cat.jrxml

sed -i 's/textAdjust="StretchHeight"/isStretchWithOverflow="true"/g' /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/LungCancer_Report_v1_cat.jrxml 
./jasperstarter compile /home/bdelolmo/Escriptori/PIPELINE_CANCER/BIN_FOLDER/JASPERREPORTS/MyReports/cat/LungCancer_Report_v1_cat.jrxml