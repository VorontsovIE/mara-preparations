PATTERN = /^>(?<name>.+)::(?<chr>chr\w+):(?<from>\d+)-(?<to>\d+)\((?<strand>[+-])\)$/
ARGF.readlines.map(&:chomp).each_slice(2){|hdr,score|
  m = hdr.match(PATTERN)
  chr, from, to, name, strand = m.named_captures.values_at('chr', 'from', 'to', 'name', 'strand')
  score = score.split("\t").first
  info = [chr, from, to, name, score, strand]
  puts(info.join("\t"))
}
