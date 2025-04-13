require 'optparse'

# 50u — 50nt upstream
# 20d — 20nt downstream
# 0   — 0nt
def parse_distance(s)
  if s == '0'
    0
  else
    match = s.downcase.match(/^(\d+)([ud])$/)
    raise "distance doesn't match format"  if !match
    (match[2] == 'u') ? -Integer(match[1]) : Integer(match[1])
  end
end

# The input file is named .bed, but instead of score it has peak center/summit, columns 7+ are arbitrary

# 1 #Chr  chr1
# 2 Start 100038029
# 3 End 100038283
# 4 ClusterID chr1_100038029_100038283_fw
# 5 Center  100038155
# 6 Strand  +
# 7 GeneID  ENSG00000156875
# 8 TotalCounts 133428
# 9 MedianTPM 17.291171741273
# 10  SoftClipScore 0.558488551432752

# or 

# 1  chr1
# 2  65418
# 3  65419
# 4  hg38_v1_chr1_+_65419_65419
# 5  65418
# 6  +
def parse_line(line, position_getter: ->(row){ Integer(row[4]) })
  row = line.split("\t")

  {
    chr: row[0], from: Integer(row[1]), to: Integer(row[2]),
    name: row[3],
    position: position_getter.call(row),
    strand: row[5],
  }
end

center_column = 5 # 1-based
option_parser = OptionParser.new{|opts|
  opts.on('--center-column NUM'){|val|
    center_column = Integer(val)
  }
}

option_parser.parse!(ARGV)

flank_5 = parse_distance(ARGV[0])
flank_3 = parse_distance(ARGV[1])

input_fn = ARGV[2]  # ./stages/stage_01/Credible_TSSclusters_regions.bed

lines = ['stdin', '-'].include?(input_fn.downcase) ? STDIN.readlines : File.readlines(input_fn)
regions = lines.reject{|l|
  l.start_with?('#')
}.map{|l|
  parse_line(l, position_getter: ->(row){ Integer(row[center_column - 1]) }) # center_column is 1-based
}

regions.map{|region|
  case region[:strand]
  when '+'
    from = region[:position] + flank_5
    to   = region[:position] + flank_3 + 1
  when '-'
    from = region[:position] - flank_3
    to   = region[:position] - flank_5 + 1
  else
    raise "Unknown strand #{strand}"
  end
  info = [region[:chr], from, to, region[:name], '.', region[:strand]]
  puts info.join("\t")
}
