require 'parallel'
require 'thread'

class Chapter1

  BASES = ['A','C','G','T']
  COMPLEMENTS = {
    'A' => 'T',
    'C' => 'G',
    'G' => 'C',
    'T' => 'A'
  }

  @frequent_base_arrays = {}
  @frequent_base_maps = {}

  def self.pattern_count(text, pattern)
    count = 0
    text.split('').each_cons(pattern.length) do |slice|
      s = slice.join()
      count += 1 if s == pattern
    end
    count
  end

  def self.slow_frequent_words(text, k)
    count = {}
    max_count = 0
    text.split('').each_cons(k) do |slice|
      s = slice.join()
      if count[s]
        count[s] = count[s] + 1
      else
        count[s] = 1
      end
      max_count = count[s] if count[s] > max_count
    end

    count.reject{|k,v| v < max_count}.sort_by{|k,v| v}.reverse.map{|k,v| k}
  end

  def self.frequency_base_array(k)
    unless @frequent_base_arrays[k]
      a_of_a = []
      (k-1).times{a_of_a << BASES}
      i = 0
      @frequent_base_maps[k] = {}
      @frequent_base_arrays[k] = BASES.product(*a_of_a).map{|p| kmer=p.join; @frequent_base_maps[k][kmer]=i; i+=1; kmer}
    end
    @frequent_base_arrays[k]
  end

  def self.pattern_to_number(pattern)
    fba = frequency_base_array(pattern.length)
    #fba.find_index(pattern)
    @frequent_base_maps[pattern.length][pattern]
  end

  def self.number_to_pattern(index, k)
    fba = frequency_base_array(k)
    fba[index]
  end

  def self.compute_patterns(text, k)
    fba = frequency_base_array(k)
    patterns = []
    for i in 0..(text.length - k)
      pattern   = text[i,k]
      index     = pattern_to_number(pattern)
      patterns << {kmer:pattern, index:index, count:1}
    end
    patterns
  end

  def self.frequent_words(text, k)
    frequent_patterns = []
    patterns = compute_patterns(text, k)

    # Sort by index
    patterns.sort_by!{|hsh| hsh[:index]}

    # Count
    for i in 1..(patterns.length-1)
      if patterns[i][:index] == patterns[i-1][:index]
        patterns[i][:count] = patterns[i-1][:count] + 1
      end
    end

    max = patterns.max_by{|hsh| hsh[:count]}[:count]

    for i in 0..(text.length - k)
      if patterns[i][:count] == max
        pattern = number_to_pattern(patterns[i][:index], k)
        frequent_patterns << pattern
      end
    end

    frequent_patterns
  end


  def self.frequent_words_mismatch(text, k, d, reverse=false)

    neighborhoods = []
    for i in 0..(text.length-k)
      pattern = text[i,k]
      neighborhoods += neighbors(pattern, d)
    end

    patterns = []
    for i in 0..(neighborhoods.length-1)
      pattern   = neighborhoods[i]
      pattern_r = reverse_complement(pattern)
      index = pattern_to_number(pattern)
      patterns << {kmer:pattern, index:index, count:1}
      patterns << {kmer:pattern_r, index:pattern_to_number(pattern_r), count:1} if reverse
    end

    # Sort by index
    patterns.sort_by!{|hsh| hsh[:index]}

    # Count
    counts = {}
    for i in 1..(patterns.length-1)
      if patterns[i][:index] == patterns[i-1][:index]
        patterns[i][:count] = patterns[i-1][:count] + 1

        counts[patterns[i][:kmer]] = patterns[i][:count]
      else
        counts[patterns[i][:kmer]] = 1
      end
    end

    # Calculate Max. Has to be greater than 1
    max = counts.max_by{|k,v| v}[1]
    max = 2 if max < 2

    frequent_patterns = []

    counts.each do |k,v|
      if v == max
        frequent_patterns << k
        puts "#{k} = #{v}"
      end
    end

    frequent_patterns
  end

  def self.skew(genome)
    skew_array = []
    skew = 0
    skew_array << skew
    for i in 0..(genome.length-1)
      nucleotide = genome[i,1]
      if nucleotide == 'C'
        skew = skew - 1
      elsif nucleotide == 'G'
        skew = skew + 1
      end
      skew_array << skew
    end
    skew_array
  end

  def self.skew_file(source, target, skip=1)
    genome = read_fasta(source)
    puts 'file read'
    skew_array = skew(genome)
    File.open(target, 'w') do |file|
      skew_array.each_with_index do |v,index|
        file.puts "#{index},#{v}" if index % skip == 0
      end
    end
    nil
  end

  def self.get_window(genome,o,l, target)
    text = genome[o,l]
    File.open(target, 'w') do |file|
      file.write text
    end
    nil
  end

  def self.read_fasta(filename)
    text = ''
    f = File.open(filename,'r')
    f.each_line do |line|
      next if line[0,1] == '>'
      text += line.chomp
    end
    f.close
    text
  end

  def self.all_min(skew_array)
    min = skew_array.min
    skew_array.each_index.select{|i| skew_array[i] == min}
  end

  def self.hamming_distance(text1, text2)
    distance = 0
    for i in 0..(text1.length-1)
      distance += 1 if text1[i,1] != text2[i,1]
    end
    distance
  end

  def self.approximate_pattern_count(text, pattern, mismatches=0)
    count = 0
    k = pattern.length
    text.split('').each_cons(k) do |slice|
      s = slice.join()
      hd = hamming_distance(s,pattern)
      count += 1 if hd <= mismatches
    end
    count
  end

  def self.neighbors(pattern, d)
    patterns = []
    patterns << pattern

    if d > 0
      chars = pattern.split('')

      d = d - 1

      for pos in 0..(chars.length-1)

        new_chars = chars.dup
        changes = BASES - [ chars[pos] ]

        changes.each do |change|
          new_chars[pos] = change
          new_pattern = new_chars.join

          patterns << new_pattern

          if d > 0
            patterns += neighbors(new_pattern, d)
          end

        end

      end
    end

    patterns.uniq
  end

  def self.reverse_complement(text)
    text.split('').reverse.map{|c| COMPLEMENTS[c] }.join
  end


end