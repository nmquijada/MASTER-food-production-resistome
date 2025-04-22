

`ls "#{ARGV[0]}" > list.txt`
m=''
m<<`wc -l list.txt
m=m.split("\s")[0]
`mkdir viromeqc`


n=1
aa=File.open("list.txt").each_line do |line|
line.chomp!
puts "...runing viromeQC analysis on sample #{line} (#{n}/#{m})"
`python3.8 viromeQC.py -i "#{ARGV[0]}"/"#{line}"/"#{line}"_R1.fastq.bz2 -o viromeqc/report_"#{line}".txt --bowtie2_threads 64 --diamond_threads 64 --minlen 0 --minqual 0`
n+=1
end
aa.close
