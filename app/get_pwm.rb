require 'optparse'

def read_matrix(fn, num_columns: 4)
  if fn == 'stdin'
    lines = $stdin.readlines
    name ||= 'motif'
  else
    lines = (fn == 'stdin') ? $stdin.readlines : File.readlines(fn)
    name = File.basename(fn, File.extname(fn))
  end
  lines = lines.map(&:strip)
  rows = lines.map{|l| l.split }
  unless (rows[0].size == num_columns) && rows[0].all?{|x| Float(x, exception: false) }
    hdr = lines.first
    rows.shift
    name = (hdr.start_with?('>') ? hdr[1..-1].strip : hdr.strip).split.first
  end
  matrix = rows.map{|row|
    row.map{|x| Float(x) }
  }
  raise  if matrix.empty?
  raise  unless matrix.all?{|row| row.size == num_columns }
  {name: name, matrix: matrix}
end

#############################

def matrix_as_string(model, transpose_output: false)
  res = [">#{model[:name]}"]
  matrix = transpose_output ? model[:matrix].transpose : model[:matrix]
  res += matrix.map{|row| row.join("\t") }
  res.join("\n")
end

#############################

def calculate_pseudocount(count, pseudocount: :log)
  case pseudocount
  when :log
    Math.log([count, 2].max);
  when :sqrt
    Math.sqrt(count)
  else Numeric
    pseudocount
  end
end

#############################

def pcm2pfm(pcm)
  pfm_matrix = pcm[:matrix].map{|row|
    norm = row.sum
    row.map{|x| x.to_f / norm }
  }
  {name: pcm[:name], matrix: pfm_matrix}
end

def pfm2pcm(pfm, word_count: 100)
  pcm_matrix = pfm[:matrix].map{|row|
    row.map{|el| el * word_count }
  }
  {name: pfm[:name], matrix: pcm_matrix}
end

def pcm2pwm(pcm, pseudocount: :log)
  pwm_matrix = pcm[:matrix].map{|row|
    count = row.sum
    row.map{|el|
      pseudocount_value = calculate_pseudocount(count, pseudocount: pseudocount)
      numerator = el + 0.25 * pseudocount_value
      denominator = 0.25 * (count + pseudocount_value)
      Math.log(numerator / denominator)
    }
  }
  {name: pcm[:name], matrix: pwm_matrix}
end

#############################

pseudocount = :log
word_count = 1000


motif_format = nil
output_format = :pwm
transpose_output = false
option_parser = OptionParser.new{|opts|
  opts.on('--transpose-output', 'Transpose output matrix'){ transpose_output = true }
  opts.on('--pfm', '--ppm', 'Force use of PFM matrix'){ motif_format = :pfm }
  opts.on('--pcm', 'Force use of PCM matrix'){ motif_format = :pcm }
  opts.on('--to-pfm', 'Convert not to PWM but to PFM') { output_format = :pfm }
  opts.on('--to-pcm', 'Convert not to PWM but to PCM') { output_format = :pcm }
  opts.on('--pseudocount VALUE', 'PCM --> PWM pseudocount'){|value| pseudocount = Float(value) rescue value.to_sym }
  opts.on('--word-count VALUE', 'PFM --> PCM word count'){|value| word_count = Float(value) }
}
option_parser.parse!(ARGV)

filename = ARGV[0]
motif_format = File.extname(filename)
if motif_format == '.pcm'
  motif_format = :pcm
elsif motif_format == '.pfm' || motif_format == '.ppm'
  motif_format = :pfm
elsif motif_format == '.pwm'
  motif_format = :pwm
end
option_parser.parse!(ARGV)

if output_format == :pwm
  if motif_format == :pfm
    pfm = read_matrix(filename, num_columns: 4)
    pcm = pfm2pcm(pfm, word_count: word_count)
    pwm = pcm2pwm(pcm, pseudocount: pseudocount)
  elsif motif_format == :pcm
    pcm = read_matrix(filename, num_columns: 4)
    pwm = pcm2pwm(pcm, pseudocount: pseudocount)
  elsif motif_format == :pwm
    pwm = read_matrix(filename, num_columns: 4)
  else
    raise "Unknown motif format `#{motif_format}`"
  end
  puts matrix_as_string(pwm, transpose_output: transpose_output)
elsif output_format == :pcm
  if motif_format == :pfm
    pfm = read_matrix(filename, num_columns: 4)
    pcm = pfm2pcm(pfm, word_count: word_count)
  elsif motif_format == :pcm
    pcm = read_matrix(filename, num_columns: 4)
    pfm = pcm2pfm(pcm)
    pcm = pfm2pcm(pfm, word_count: word_count) # renormalize
  elsif motif_format == :pwm
    raise "Can't convert PWM to PCM"
  else
    raise "Unknown motif format `#{motif_format}`"
  end
  puts matrix_as_string(pcm, transpose_output: transpose_output)
elsif output_format == :pfm
  if motif_format == :pfm
    pfm = read_matrix(filename, num_columns: 4)
  elsif motif_format == :pcm
    pcm = read_matrix(filename, num_columns: 4)
    pfm = pcm2pfm(pcm)
  elsif motif_format == :pwm
    raise "Can't convert PWM to PFM"
  else
    raise "Unknown motif format `#{motif_format}`"
  end
  puts matrix_as_string(pfm, transpose_output: transpose_output)
else
  raise "Unknown output format `#{output_format}`"
end

