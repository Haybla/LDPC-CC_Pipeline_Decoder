# LDPC-CC_Pipeline_Decoder

Copyright (c) 2014-2015 Mokyy and Haybla. All rights reserved.

This is the source codes of a CUDA application named LDPC-CC_Pipeline_Decoder. 
Techniques of design and optimization can be found in the article "High Throughput 
Pipeline Decoder for LDPC Convolutional Codes on GPU" published in the IEEE 
Communications Letters journal (The article was accepted just now).

In order to compile this project, Linux OS is recommended. Just open a terminal 
and go into the "src" directory. 

>cd src

Then compile all the codes using "make".

>make

An executed file named "decode" will be generated if everything goes right. Just 
execute the "decode" as follows.

>./decode

And you will see some information about decoding times and throughputs printed in 
the terminal. You can clean all compilation using "make clean".

>make clean

#Besides
Macro definitions are written in the file named "totalDefine.h". You can enable/disable 
some of them to change the mode of the decoder. Descriptions are as follows.

CODE1: represents the LDPC-CC constructed by (4096, 10240) LDPC codes in CCSDS standard. 

CODE2: represents the LDPC-CC constructed by (4096, 7168) LDPC codes in CCSDS standard.

LINUX: once defined, the project must be compiled and run on Linux OS.

TEST_PERF： once defined, the project can simulate BER&FER performance of the LDPC-CC.

