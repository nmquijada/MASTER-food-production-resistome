
if ARGV[0]==nil
puts ""
puts ".........................................................................."
puts "... Sorry! you forgot to add the folder name which contains .bam files ..."
puts "...          try again adding this information, please                 ..."
puts ".........................................................................."
puts ""
exit
elsif ARGV[0]=~/\//
`ls "#{ARGV[0]}"* > list.txt`
else 
`ls "#{ARGV[0]}"/* > list.txt`
end

l=''
l << `wc -l list.txt`
l=l.split("\s")[0]
m=0

sample=''
out2=File.new("report_different_match.txt","w")
out3=File.new("report_different_family.txt","w")
problems=File.new("problem5gene_#{ARGV[1]}.txt","w")

hfam={}
genes=[]
hgenes={}
hgen={}
hfamgen={}
gen=''
n=0
gg=File.open("gene_list.txt").each_line do |line|
line.chomp!
#Family  Genefam Gene    Hit
#Aminoglycoside  ant(2'')-Ia     ant(2'')-Ia     ant(2'')-Ia_1_X04555
if n!=0
col=line.split("\t")
gen=col[3].gsub("\s","") #to avoid space problems
hfam[gen]=col[0]
hfamgen[gen]=col[1]
hgen[gen]=col[2]
genes << gen
end
n+=1
end
gg.close

samples=[]
hsamples={}
ll=File.open("list.txt").each_line do |line|
line.chomp!
if line =~ /\/(.*)\.\w\w\w/
samples << $1
end
end
ll.close

genes.each_index{|a| hgenes[genes[a]]=a}
samples.each_index{|b| hsamples[samples[b]]=b}
matrix=Array.new(genes.length()){|i| Array.new(samples.length()) {|x| 0}}

#puts samples
gene=''
aa=File.open("list.txt").each_line do |file|
file.chomp!
if file =~ /\/(.*).bam/ or file =~ /\/(.*).sam/ 
sample=$1
m+=1
puts "...processing sample #{sample} (#{m}/#{l})"
`samtools sort -n "#{file}" > file1.bam`
`samtools view file1.bam > file2.sam`
read=''
        bb=File.open("file2.sam").each_line do |line|
 	line.chomp!
	cols=line.split("\t")
gene=cols[2]

	if cols[0]!=read
puts sample
#puts cols[0]
puts gene

	matrix[hgenes[gene]][hsamples[sample]]+=1
	read=cols[0]
        elsif cols[6]!="="
	out2.puts line
		if hfam[gene]!=hfam[cols[6]]
		out3.puts "#{hfam[gene]}\t#{hfam[cols[6]]}\t#{line}"
puts gene
puts hgenes[gene]
puts sample
puts hsamples[sample]
		matrix[hgenes[gene]][hsamples[sample]]+=1		
		end
	end
   end
   bb.close
end
end
aa.close

name=''
#if ARGV[0] =~ /sam_files_(\S+)/
#name=$1.gsub("/","")
matfile=File.new("matrix_ResFinder.txt","w")
#end

mm=-1
matfile.print "Antibiotic_Family\tGenefam\tGene\tHit"
samples.each {|i| matfile.print "\t#{i}"}
matfile.print "\n"
matrix.each {|i| mm+=1
matfile.print "#{hfam[genes[mm]]}\t#{hfamgen[genes[mm]]}\t#{hgen[genes[mm]]}\t#{genes[mm]}"
i.each {|x| matfile.print "\t#{x}"}
matfile.print "\n"}
