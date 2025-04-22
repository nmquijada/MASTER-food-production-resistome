
# The reads has to been in a folder named "reads" in hte same pathway as the script, but you can change pathway or name on line 7

`mkdir sam_files`
`ls "#{ARGV[0]}" > list.txt`


aa=File.open("list.txt").each_line do |line|
line.chomp!
puts "...running resfinder search on sample #{line}"
`bowtie2 -p 64 --no-unal --very-sensitive --end-to-end -x bt2_resfinderdb -1 "#{ARGV[0]}"/"#{line}"/"#{line}"_R1.fastq.bz2 -2 "#{ARGV[0]}"/"#{line}"/"#{line}"_R2.fastq.bz2 -S sam_files/"#{line}".sam`
end
aa.close

