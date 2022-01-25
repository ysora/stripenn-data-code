stripenn compute --cool /mnt/data0/sora/stripe/hic/Arima-HiC_HPAP-018_Normal_PLN_CD4_1M_ultraDeep_200885694_S3_001_genome.bwt2pairs.validPairs.mcool::resolutions/5000 \
--norm VC_SQRT \
-n 10 \
--out /mnt/data0/sora/stripe/stripenn_out/human_mouse_conserve/human_CD4_test \
-w 16 \
--canny 2.0 \
-k 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X \
--maxpixel 0.97,0.98,0.99
