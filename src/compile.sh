which perl || { echo "perl is not found"; exit 1; }
echo "perl is ok"
which python || { echo "python is not found"; exit 1; }
echo "python is ok"
echo "installing bowtie2"
if [ ! -d "../bin/bowtie2-2.2.1/" ]; then
  echo "Fold ../bin/bowtie2-2.2.1/ does not exist."
  exit 1
fi
cd ../bin/bowtie2-2.2.1/ || { exit 1; }
make || { echo "make bowtie2 failed"; exit 1; }
echo "Done"
