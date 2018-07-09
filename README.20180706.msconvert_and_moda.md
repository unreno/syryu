#	MSConvert and MODa

##	20180706


1. Connect to ftp://massive.ucsd.edu/MSV000079053. There are over 100 bacteria in this folder. 
2. Download fasta files and raw files for each bacteria. 
    a. fasta files: Under sequence folder, there are many folders (one folder for one bacteria). In bacteria folder, download files ending with .fasta.
    b. RAW files: Under raw folder, there are many folders (one folder for one bacteria). In bacteria folder, download files ending with RAW.  
3. Convert RAW to mgf.
    a. Download ProteoWizard. http://proteowizard.sourceforge.net/download.html
        May have to download window version to convert raw files. 
    b. Convert RAW to mgf using msconvert command
    http://proteowizard.sourceforge.net/tools/msconvert.html
    ex) msconvert data.RAW --mgf --filter "msLevel 2" --filter "zeroSample removeExtra"
4. Run MODa. I attached an example parameter file (foo.txt) and a tutorial (README.pdf). For each mgf file, please change "Spectra=" and "Fasta=" and run MODa. 
5. Transfer results files to Azure storage. 

Please let me know if you have questions!
Thank you, Jake.





Seems msconvert.exe needs MSFileReader to read RAW files.

Still not working? Adding XcalDLL

convertion seems to be working.





C:\Users\jake>\Ruby25-x64\bin\ruby.exe msconvert.rb




START /B \Ruby25-x64\bin\ruby.exe msconvert.rb













------




for f in Acidiphilium_cryptum_JF-5 Actinosynnema_mirum_DSM_43827 Anabaena_variabilis Anaeromyxobacter_dehalogenans Anaplasma_phagocytophilium Arthrobacter_sp_FB24 Bacillus_anthracis_Ames Bacillus_anthracis_Sterne Bacillus_subtilis_168 Bartonella_henselae_Houston-1 Borrelia_burdorferi_B31 Brachybacterium_faecium_DSM_4810 Burkholderia_mallei Candidatus_chloracidobacterium_thermophilum Caulobacter_crescentus_CB15 Cellulomonas_flavigena_DSM_20109 Cenarchaeum_symbiosum Chlorobaculum_tepidum_WT Chloroflexus_aurantiacus Clostridium_thermocellum Cryptobacterium_curtum_DSM_15641 Cyanobacterium_synechocystis_PCC6803 Cyanothece_sp_ATCC51142 Cyanothece_strain_ATCC51472 Cyanothece_strain_PCC7424 Cyanothece_strain_PCC7425 Cyanothece_strain_PCC7822 Cyanothece_strain_PCC8801 Cyanothece_strain_PCC8802 Dehalococcoides_ethenogenes Deinococcus_radiodurans_R1 Delta_proteobacterium_NaphS2 Desulfovibrio_desulfuricans_G20 Desulfovibrio_sp_ND132 Desulfovibrio_vulgaris_Hildenborough Dethiosulfovibrio_peptidovorans_DSM_11002 Ehrlichia_chaffeensis Enterobacter_cloacae_SCF1 Escherichia_coli_BL21 Escherichia_coli_K-12 Escherichia_coli_RK4353 Fibrobacter_succinogenes_S85 Geobacter_bemidjiensis_Bem_T Geobacter_metallireducens_GS-15 Geobacter_sulfurreducens_PCA Geobacter_uraniumreducens Haloferax_volcanii Halogeometricum_borinquense_DSM_11551 Halorhabdus_utahensis_DSM_12940 Heliobacterium_modesticaldum Kineococcus_radiotolerans_SRS30216 Kosmotoga_olearia_TBF_19-5-1 Magnetospirillum_magneticum_AMB-1 Methanosarcina_barkeri Methanospirillum_hungatei_JF-1 Methylophilales_HTCC2181 Mycobacterium_tuberculosis Nakamurella_multipartita_DSM_44233 Nocardiopsis_dassonvillei_DSM_43111 Novosphingobium_aromaticivorans_F199 Opitutaceae_bacterium_TAV2 Pelagibacter_ubique_HTC1062 Pelobacter_carbinolicus_DSM_2380 Prochlorococcus Pseudomonas_aerunginosa Pseudomonas_fluorescens_PfO-1 Pseudonocardia_sp Ralstonia_pickettii Rhodobacter_capsulatus_SB1003 Rhodobacter_sphaeroides_2.4.1 Rhodopseudomonas_palustris Roseiflexus_castenholzii Saccharomonospora_viridis_DSM_43017 Salmonella_typhi_TY2 Salmonella_typhimurium_ATCC_14028 Salmonella_typhimurium_LT2 Sanguibacter_keddieii_DSM_10542 Shewanella_amazonensis_SB2B Shewanella_baltica_OS155 Shewanella_baltica_OS185 Shewanella_baltica_OS195 Shewanella_baltica_OS223 Shewanella_denitrificans_OS217 Shewanella_frigidimarina_NCIMB_400 Shewanella_loihica_PV-4 Shewanella_oneidensis_MR-1 Shewanella_putrefaciens_200 Shewanella_putrefaciens_CN-32 Shewanella_putrefaciens_W3-18-1 Shewanella_sp_ANA-3 Shewanella_sp_MR-4 Shewanella_sp_MR-7 Sinorhizobium_medicae Sinorhizobium_meliloti_1021 Slackia_heliotrinireducens_DSM_20476 Stackebrandtia_nassauensis_DSM_44728 Sulfolobus_acidocaldarius_DSM_639 Synechococcus_sp_PCC7002 Syntrophobacter_fumaroxidans Thermobispora_bispora_DSM_43833 Thermosynechococcus_elongatus_BP-1 Thermosynechococcus_sp_NAK55 Thermotoga_maritima Thiocapsa_marina_DSM_5653T Verrucomicrobium_sp_TAV1 Verrucomicrobium_sp_TAV5 Xylanimonas_cellulosilytica_DSM_15894 Yersinia_enterocolitica Yersinia_pestis_CO92 Yersinia_pestis_KIM Yersinia_pestis_Pestoides_F Yersinia_pseudotuberculosis_IP_32953 Yersinia_pseudotuberculosis_PB1_Plus ; do

echo "START /B \Ruby25-x64\bin\ruby.exe msconvert.rb $f > $f.log"
done



