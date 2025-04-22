cluster=''
group=''
hhit={}
hgen={}
hfam={}
shortname=''

aa=File.open(ARGV[0]).each_line do |line|
line.chomp!
#ant(2'')-Ia_1_X04555   Aminoglycoside
col=line.split("\t")
if line =~ /\w/ # to avoid problems with empty lines
if col[0].length > 19 # it is because CDHIT put maximum 19 characters on .clstr output file, and we know that these are different for all genes/hits
hhit[col[0][0,19]]=col[0]
hgen[col[0][0,19]]=col[0].split("\_")[0]
hfam[col[0][0,19]]=col[1]
else
hhit[col[0]]=col[0]
hgen[col[0]]=col[0].split("\_")[0]
hfam[col[0]]=col[1]
#puts col[0].length
end
end # of \w line
end
aa.close

puts hgen

out=File.new("genes_90i.txt","w")
out.puts "Cluster\tAntibiotic_family\tGene\tHit"

bb=File.open(ARGV[1]).each_line do |line|
line.chomp!
if line =~ /^>(\S+)\s+(\S+)/
cluster="#{$1}#{$2}"
#elsif line =~ />(\S+)\_\d+\_/ # and line =~ /\*/
elsif line =~ />(\S+)\.\.\./ # and line =~ /\*/
group=$1
out.puts "#{cluster}\t#{hfam[$1]}\t#{hgen[$1]}\t#{hhit[$1]}" #\t#{$1}"
end
end
bb.close
