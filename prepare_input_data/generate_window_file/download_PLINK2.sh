cd /software/team152/jingling_software
rm -rf plink*
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip
unzip -o plink2_linux_avx2_20250129.zip
ls -lh /software/team152/jingling_software
unzip -l plink2_linux_avx2_20250129.zip
export PATH="/software/team152/jingling_software:$PATH"
chmod +x plink2
