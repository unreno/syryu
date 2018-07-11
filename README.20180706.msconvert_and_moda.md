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

echo "START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb $f > $f.log 2>&1"

done






START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Acidiphilium_cryptum_JF-5 > Acidiphilium_cryptum_JF-5.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Actinosynnema_mirum_DSM_43827 > Actinosynnema_mirum_DSM_43827.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Anabaena_variabilis > Anabaena_variabilis.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Anaeromyxobacter_dehalogenans > Anaeromyxobacter_dehalogenans.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Anaplasma_phagocytophilium > Anaplasma_phagocytophilium.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Arthrobacter_sp_FB24 > Arthrobacter_sp_FB24.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Bacillus_anthracis_Ames > Bacillus_anthracis_Ames.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Bacillus_anthracis_Sterne > Bacillus_anthracis_Sterne.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Bacillus_subtilis_168 > Bacillus_subtilis_168.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Bartonella_henselae_Houston-1 > Bartonella_henselae_Houston-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Borrelia_burdorferi_B31 > Borrelia_burdorferi_B31.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Brachybacterium_faecium_DSM_4810 > Brachybacterium_faecium_DSM_4810.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Burkholderia_mallei > Burkholderia_mallei.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Candidatus_chloracidobacterium_thermophilum > Candidatus_chloracidobacterium_thermophilum.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Caulobacter_crescentus_CB15 > Caulobacter_crescentus_CB15.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cellulomonas_flavigena_DSM_20109 > Cellulomonas_flavigena_DSM_20109.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cenarchaeum_symbiosum > Cenarchaeum_symbiosum.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Chlorobaculum_tepidum_WT > Chlorobaculum_tepidum_WT.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Chloroflexus_aurantiacus > Chloroflexus_aurantiacus.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Clostridium_thermocellum > Clostridium_thermocellum.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cryptobacterium_curtum_DSM_15641 > Cryptobacterium_curtum_DSM_15641.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanobacterium_synechocystis_PCC6803 > Cyanobacterium_synechocystis_PCC6803.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_sp_ATCC51142 > Cyanothece_sp_ATCC51142.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_strain_ATCC51472 > Cyanothece_strain_ATCC51472.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_strain_PCC7424 > Cyanothece_strain_PCC7424.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_strain_PCC7425 > Cyanothece_strain_PCC7425.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_strain_PCC7822 > Cyanothece_strain_PCC7822.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_strain_PCC8801 > Cyanothece_strain_PCC8801.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Cyanothece_strain_PCC8802 > Cyanothece_strain_PCC8802.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Dehalococcoides_ethenogenes > Dehalococcoides_ethenogenes.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Deinococcus_radiodurans_R1 > Deinococcus_radiodurans_R1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Delta_proteobacterium_NaphS2 > Delta_proteobacterium_NaphS2.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Desulfovibrio_desulfuricans_G20 > Desulfovibrio_desulfuricans_G20.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Desulfovibrio_sp_ND132 > Desulfovibrio_sp_ND132.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Desulfovibrio_vulgaris_Hildenborough > Desulfovibrio_vulgaris_Hildenborough.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Dethiosulfovibrio_peptidovorans_DSM_11002 > Dethiosulfovibrio_peptidovorans_DSM_11002.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Ehrlichia_chaffeensis > Ehrlichia_chaffeensis.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Enterobacter_cloacae_SCF1 > Enterobacter_cloacae_SCF1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Escherichia_coli_BL21 > Escherichia_coli_BL21.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Escherichia_coli_K-12 > Escherichia_coli_K-12.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Escherichia_coli_RK4353 > Escherichia_coli_RK4353.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Fibrobacter_succinogenes_S85 > Fibrobacter_succinogenes_S85.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Geobacter_bemidjiensis_Bem_T > Geobacter_bemidjiensis_Bem_T.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Geobacter_metallireducens_GS-15 > Geobacter_metallireducens_GS-15.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Geobacter_sulfurreducens_PCA > Geobacter_sulfurreducens_PCA.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Geobacter_uraniumreducens > Geobacter_uraniumreducens.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Haloferax_volcanii > Haloferax_volcanii.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Halogeometricum_borinquense_DSM_11551 > Halogeometricum_borinquense_DSM_11551.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Halorhabdus_utahensis_DSM_12940 > Halorhabdus_utahensis_DSM_12940.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Heliobacterium_modesticaldum > Heliobacterium_modesticaldum.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Kineococcus_radiotolerans_SRS30216 > Kineococcus_radiotolerans_SRS30216.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Kosmotoga_olearia_TBF_19-5-1 > Kosmotoga_olearia_TBF_19-5-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Magnetospirillum_magneticum_AMB-1 > Magnetospirillum_magneticum_AMB-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Methanosarcina_barkeri > Methanosarcina_barkeri.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Methanospirillum_hungatei_JF-1 > Methanospirillum_hungatei_JF-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Methylophilales_HTCC2181 > Methylophilales_HTCC2181.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Mycobacterium_tuberculosis > Mycobacterium_tuberculosis.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Nakamurella_multipartita_DSM_44233 > Nakamurella_multipartita_DSM_44233.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Nocardiopsis_dassonvillei_DSM_43111 > Nocardiopsis_dassonvillei_DSM_43111.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Novosphingobium_aromaticivorans_F199 > Novosphingobium_aromaticivorans_F199.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Opitutaceae_bacterium_TAV2 > Opitutaceae_bacterium_TAV2.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Pelagibacter_ubique_HTC1062 > Pelagibacter_ubique_HTC1062.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Pelobacter_carbinolicus_DSM_2380 > Pelobacter_carbinolicus_DSM_2380.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Prochlorococcus > Prochlorococcus.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Pseudomonas_aerunginosa > Pseudomonas_aerunginosa.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Pseudomonas_fluorescens_PfO-1 > Pseudomonas_fluorescens_PfO-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Pseudonocardia_sp > Pseudonocardia_sp.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Ralstonia_pickettii > Ralstonia_pickettii.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Rhodobacter_capsulatus_SB1003 > Rhodobacter_capsulatus_SB1003.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Rhodobacter_sphaeroides_2.4.1 > Rhodobacter_sphaeroides_2.4.1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Rhodopseudomonas_palustris > Rhodopseudomonas_palustris.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Roseiflexus_castenholzii > Roseiflexus_castenholzii.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Saccharomonospora_viridis_DSM_43017 > Saccharomonospora_viridis_DSM_43017.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Salmonella_typhi_TY2 > Salmonella_typhi_TY2.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Salmonella_typhimurium_ATCC_14028 > Salmonella_typhimurium_ATCC_14028.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Salmonella_typhimurium_LT2 > Salmonella_typhimurium_LT2.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Sanguibacter_keddieii_DSM_10542 > Sanguibacter_keddieii_DSM_10542.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_amazonensis_SB2B > Shewanella_amazonensis_SB2B.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_baltica_OS155 > Shewanella_baltica_OS155.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_baltica_OS185 > Shewanella_baltica_OS185.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_baltica_OS195 > Shewanella_baltica_OS195.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_baltica_OS223 > Shewanella_baltica_OS223.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_denitrificans_OS217 > Shewanella_denitrificans_OS217.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_frigidimarina_NCIMB_400 > Shewanella_frigidimarina_NCIMB_400.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_loihica_PV-4 > Shewanella_loihica_PV-4.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_oneidensis_MR-1 > Shewanella_oneidensis_MR-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_putrefaciens_200 > Shewanella_putrefaciens_200.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_putrefaciens_CN-32 > Shewanella_putrefaciens_CN-32.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_putrefaciens_W3-18-1 > Shewanella_putrefaciens_W3-18-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_sp_ANA-3 > Shewanella_sp_ANA-3.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_sp_MR-4 > Shewanella_sp_MR-4.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Shewanella_sp_MR-7 > Shewanella_sp_MR-7.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Sinorhizobium_medicae > Sinorhizobium_medicae.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Sinorhizobium_meliloti_1021 > Sinorhizobium_meliloti_1021.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Slackia_heliotrinireducens_DSM_20476 > Slackia_heliotrinireducens_DSM_20476.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Stackebrandtia_nassauensis_DSM_44728 > Stackebrandtia_nassauensis_DSM_44728.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Sulfolobus_acidocaldarius_DSM_639 > Sulfolobus_acidocaldarius_DSM_639.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Synechococcus_sp_PCC7002 > Synechococcus_sp_PCC7002.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Syntrophobacter_fumaroxidans > Syntrophobacter_fumaroxidans.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Thermobispora_bispora_DSM_43833 > Thermobispora_bispora_DSM_43833.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Thermosynechococcus_elongatus_BP-1 > Thermosynechococcus_elongatus_BP-1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Thermosynechococcus_sp_NAK55 > Thermosynechococcus_sp_NAK55.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Thermotoga_maritima > Thermotoga_maritima.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Thiocapsa_marina_DSM_5653T > Thiocapsa_marina_DSM_5653T.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Verrucomicrobium_sp_TAV1 > Verrucomicrobium_sp_TAV1.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Verrucomicrobium_sp_TAV5 > Verrucomicrobium_sp_TAV5.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Xylanimonas_cellulosilytica_DSM_15894 > Xylanimonas_cellulosilytica_DSM_15894.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Yersinia_enterocolitica > Yersinia_enterocolitica.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Yersinia_pestis_CO92 > Yersinia_pestis_CO92.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Yersinia_pestis_KIM > Yersinia_pestis_KIM.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Yersinia_pestis_Pestoides_F > Yersinia_pestis_Pestoides_F.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Yersinia_pseudotuberculosis_IP_32953 > Yersinia_pseudotuberculosis_IP_32953.log 2>&1
START /B \Ruby25-x64\bin\ruby.exe \ryulab\syryu\archive\msconvert_moda.rb Yersinia_pseudotuberculosis_PB1_Plus > Yersinia_pseudotuberculosis_PB1_Plus.log 2>&1




3.8G    Acidiphilium_cryptum_JF-5                                                         
6.2G    Actinosynnema_mirum_DSM_43827                                                         
26G     Anabaena_variabilis                                                         
5.5G    Anaeromyxobacter_dehalogenans                                                         
11G     Anaplasma_phagocytophilium                                                         
36G     Arthrobacter_sp_FB24                                                         
2.3G    Bacillus_anthracis_Ames                                                         
7.7G    Bacillus_anthracis_Sterne                                                         
133G    Bacillus_subtilis_168                                                         
67G     Bartonella_henselae_Houston-1                                                         
6.5G    Borrelia_burdorferi_B31                                                         
4.9G    Brachybacterium_faecium_DSM_4810                                                         
27G     Burkholderia_mallei                                                         
5.5G    Candidatus_chloracidobacterium_thermophilum                                                         
163G    Caulobacter_crescentus_CB15                                                         
6.5G    Cellulomonas_flavigena_DSM_20109                                                         
34G     Cenarchaeum_symbiosum                                                         
5.2G    Chlorobaculum_tepidum_WT                                                         
113G    Chloroflexus_aurantiacus                                                         
170G    Clostridium_thermocellum                                                         
2.8G    Cryptobacterium_curtum_DSM_15641                                                         
96G     Cyanobacterium_synechocystis_PCC6803                                                         
189G    Cyanothece_sp_ATCC51142                                                         
15G     Cyanothece_strain_ATCC51472                                                         
13G     Cyanothece_strain_PCC7424                                                         
18G     Cyanothece_strain_PCC7425                                                         
58G     Cyanothece_strain_PCC7822                                                         
15G     Cyanothece_strain_PCC8801                                                         
12G     Cyanothece_strain_PCC8802                                                         
188G    Dehalococcoides_ethenogenes                                                         
95G     Deinococcus_radiodurans_R1                                                         
8.6G    Delta_proteobacterium_NaphS2                                                         
99G     Desulfovibrio_desulfuricans_G20                                                         
9.6G    Desulfovibrio_sp_ND132                                                         
27G     Desulfovibrio_vulgaris_Hildenborough                                                         
5.1G    Dethiosulfovibrio_peptidovorans_DSM_11002                                                         
12G     Ehrlichia_chaffeensis                                                         
47G     Enterobacter_cloacae_SCF1                                                         
32G     Escherichia_coli_BL21                                                         
1.2T    Escherichia_coli_K-12                                                         
2.6G    Escherichia_coli_RK4353                                                         
60G     Fibrobacter_succinogenes_S85                                                         
328G    Geobacter_bemidjiensis_Bem_T                                                         
24G     Geobacter_metallireducens_GS-15                                                         
308G    Geobacter_sulfurreducens_PCA                                                         
22G     Geobacter_uraniumreducens                                                         
18G     Haloferax_volcanii                                                         
3.4G    Halogeometricum_borinquense_DSM_11551                                                         
4.0G    Halorhabdus_utahensis_DSM_12940                                                         
4.0G    Heliobacterium_modesticaldum                                                         
38G     Kineococcus_radiotolerans_SRS30216                                                         
16G     Kosmotoga_olearia_TBF_19-5-1                                                         
318M    Magnetospirillum_magneticum_AMB-1                                                         
4.3G    Methanosarcina_barkeri                                                         
5.7G    Methanospirillum_hungatei_JF-1                                                         
14G     Methylophilales_HTCC2181                                                         
290G    Mycobacterium_tuberculosis                                                         
3.8G    Nakamurella_multipartita_DSM_44233                                                         
5.7G    Nocardiopsis_dassonvillei_DSM_43111                                                         
16G     Novosphingobium_aromaticivorans_F199                                                         
21G     Opitutaceae_bacterium_TAV2                                                         
168G    Pelagibacter_ubique_HTC1062                                                         
2.0G    Pelobacter_carbinolicus_DSM_2380                                                         
41G     Prochlorococcus                                                         
16G     Pseudomonas_aerunginosa                                                         
28G     Pseudomonas_fluorescens_PfO-1                                                         
2.4G    Pseudonocardia_sp                                                         
9.3G    Ralstonia_pickettii                                                         
104G    Rhodobacter_capsulatus_SB1003                                                         
428G    Rhodobacter_sphaeroides_2.4.1                                                         
27G     Rhodopseudomonas_palustris                                                         
4.0G    Roseiflexus_castenholzii                                                         
4.4G    Saccharomonospora_viridis_DSM_43017                                                         
66G     Salmonella_typhi_TY2                                                         
1.8T    Salmonella_typhimurium_ATCC_14028                                                         
98G     Salmonella_typhimurium_LT2                                                         
3.8G    Sanguibacter_keddieii_DSM_10542                                                         
38G     Shewanella_amazonensis_SB2B                                                         
87G     Shewanella_baltica_OS155                                                         
14G     Shewanella_baltica_OS185                                                         
7.1G    Shewanella_baltica_OS195                                                         
4.1G    Shewanella_baltica_OS223                                                         
7.5G    Shewanella_denitrificans_OS217                                                         
6.9G    Shewanella_frigidimarina_NCIMB_400                                                         
7.2G    Shewanella_loihica_PV-4                                                         
875G    Shewanella_oneidensis_MR-1                                                         
7.4G    Shewanella_putrefaciens_200                                                         
10G     Shewanella_putrefaciens_CN-32                                                         
11G     Shewanella_putrefaciens_W3-18-1                                                         
6.8G    Shewanella_sp_ANA-3                                                         
7.4G    Shewanella_sp_MR-4                                                         
6.9G    Shewanella_sp_MR-7                                                         
22G     Sinorhizobium_medicae                                                         
80G     Sinorhizobium_meliloti_1021                                                         
4.9G    Slackia_heliotrinireducens_DSM_20476                                                         
5.8G    Stackebrandtia_nassauensis_DSM_44728                                                         
41G     Sulfolobus_acidocaldarius_DSM_639                                                         
297G    Synechococcus_sp_PCC7002                                                         
5.8G    Syntrophobacter_fumaroxidans                                                         
4.5G    Thermobispora_bispora_DSM_43833                                                         
13G     Thermosynechococcus_elongatus_BP-1                                                         
43G     Thermosynechococcus_sp_NAK55                                                         
66G     Thermotoga_maritima                                                         
8.4G    Thiocapsa_marina_DSM_5653T                                                         
3.1G    Verrucomicrobium_sp_TAV1                                                         
3.2G    Verrucomicrobium_sp_TAV5                                                         
13G     Xylanimonas_cellulosilytica_DSM_15894                                                         
15G     Yersinia_enterocolitica                                                         
100G    Yersinia_pestis_CO92                                                         
10G     Yersinia_pestis_KIM                                                         
48G     Yersinia_pestis_Pestoides_F                                                         
42G     Yersinia_pseudotuberculosis_IP_32953                                                         
54G     Yersinia_pseudotuberculosis_PB1_Plus  


