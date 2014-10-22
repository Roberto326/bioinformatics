class Chapter1

  BASES = ['A','C','G','T']

  @frequent_base_arrays = {}

  def self.pattern_count(text, pattern)
    count = 0
    text.split('').each_cons(pattern.length) do |slice|
      s = slice.join()
      count += 1 if s == pattern
    end
    count
  end

  def self.frequent_words(text, k)
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
      (k-1).times{ a_of_a << BASES}
      @frequent_base_arrays[k] = BASES.product(*a_of_a).map{|p| p.join}
    end
    @frequent_base_arrays[k]
  end

  def self.pattern_to_number(pattern)
    fba = frequency_base_array(pattern.length)
    fba.find_index(pattern)
  end

  def self.number_to_pattern(index, k)
    fba = frequency_base_array(k)
    fba[index]
  end

  def self.compute_frequencies(text, k)
    fba = frequency_base_array(k)
    fa = Array.new(fba.length){0}
    for i in 0..(text.length - k)
      slice = text[i,k]
      index = pattern_to_number(slice)
      fa[index] += 1
    end
    fa
  end

  def self.faster_frequent_words(text, k)
    fa = compute_frequencies(text, k)
    max = fa.max
    max_indices = []
    fa.each_index {|index| max_indices << index if fa[index] == max }
    max_indices.map{|index| number_to_pattern(index, k)}
  end

end