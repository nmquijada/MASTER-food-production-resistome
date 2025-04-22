
cluster=''
group=''

aa=File.open("cdhit.txt.clstr").each_line do |line|
line.chomp!
if line =~ /^>(\S+)\s+(\S+)/
cluster="#{$1}#{$2}"
#elsif line =~ />(\S+)\_\d+\_/ # and line =~ /\*/
elsif line =~ />(\S+)\.\.\./ # and line =~ /\*/
group=$1
puts "#{cluster}\t#{$1}" #\t#{$1}"
end
end
aa.close
