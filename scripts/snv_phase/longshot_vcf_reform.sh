sex=${2}
file=${1}

if [ $sex == "male" ];then
	cat ${file} | whatshap unphase - | awk -v OFS='\t' '{
		if(substr($1,1,1)=="#")
			print$0; 
		else {
			if(($1=="chrX" && $2>=2394411 && $2<=153925834) || $1=="chrY" || $1=="chrM" ) { 
				split($10,info,":");
				if(info[1]=="1/0") 
					print $1,$2,$3,$4,$5,$6,$7,$8,$9,"0/0:34:0:0/0:33.01";
				else print$0 
			}
			else print$0 
		}
	}' | awk -v OFS='\t' '{
                if(substr($1,1,1)=="#")
                        print$0;
                else {
                        split($10,info,":");
                        if(info[3]==0)
                                print $1,$2,$3,$4,$5,$6,$7,$8,$9,"./.:34:0:0/0:33.01";
                        else print$0
                }
        }'
fi


if [ $sex == "female" ];then
        cat ${file} | whatshap unphase - | awk -v OFS='\t' '{
                if(substr($1,1,1)=="#")
                        print$0;
                else {
                        if($1=="chrM") {
                                split($10,info,":");
                                if(info[1]=="1/0")
                                        print $1,$2,$3,$4,$5,$6,$7,$8,$9,"0/0:34:0:0/0:33.01";
                                else print$0
                        }
			else if($1=="chrY") next
                        else print$0
                }
        }' | awk -v OFS='\t' '{
		if(substr($1,1,1)=="#")
                        print$0;
		else {
			split($10,info,":");
			if(info[3]==0)
				print $1,$2,$3,$4,$5,$6,$7,$8,$9,"./.:34:0:0/0:33.01";
			else print$0
		}
	}'
fi

