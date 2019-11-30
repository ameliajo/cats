#/bin/sh

ls input_$1
cp input_$1 ~/NJOY2016/bin
cd ~/NJOY2016/bin
make && ./njoy < input_$1
cd -
cp ~/NJOY2016/bin/njoy_output_.py ./njoy_output_$1.py
