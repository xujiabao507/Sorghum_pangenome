#export LD_LIBRARY_PATH="/ifs4/BC_RD/USER/helijuan/software/Programm/boost_1_57_0_v2/boost_1_57_0/v2/lib/:/ifs4/BC_PUB/biosoft/pipeline/Package/GCC/gcc_4.9.0/lib64:$LD_LIBRARY_PATH"
q1="/zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S04.Platanus/step01.input/Sorghum.assembled.fastq"
q2="/zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S04.Platanus/step01.input/Sorghum.unassembled.forward.fastq /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S04.Platanus/step01.input/Sorghum.unassembled.reverse.fastq"
q3="/zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S04.Platanus/step01.input/wHAXPI051963-113_1.fastq /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S04.Platanus/step01.input/wHAXPI051963-113_2.fastq"
/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/Platanus_v1.2.4/platanus  assemble -o Sorghum -k 65 -c 3 -t 80 -m 500 -u 0.3 -f $q1 $q2 $q3
