require 'chapter1'

class Chapter2

  TYROCIDINE_B1 = 'VKLFPWFNQY'

  CODON_PROTEIN = {
    'AAA' => 'K',
    'AAC' => 'N',
    'AAG' => 'K',
    'AAU' => 'N',
    'ACA' => 'T',
    'ACC' => 'T',
    'ACG' => 'T',
    'ACU' => 'T',
    'AGA' => 'R',
    'AGC' => 'S',
    'AGG' => 'R',
    'AGU' => 'S',
    'AUA' => 'I',
    'AUC' => 'I',
    'AUG' => 'M',
    'AUU' => 'I',
    'CAA' => 'Q',
    'CAC' => 'H',
    'CAG' => 'Q',
    'CAU' => 'H',
    'CCA' => 'P',
    'CCC' => 'P',
    'CCG' => 'P',
    'CCU' => 'P',
    'CGA' => 'R',
    'CGC' => 'R',
    'CGG' => 'R',
    'CGU' => 'R',
    'CUA' => 'L',
    'CUC' => 'L',
    'CUG' => 'L',
    'CUU' => 'L',
    'GAA' => 'E',
    'GAC' => 'D',
    'GAG' => 'E',
    'GAU' => 'D',
    'GCA' => 'A',
    'GCC' => 'A',
    'GCG' => 'A',
    'GCU' => 'A',
    'GGA' => 'G',
    'GGC' => 'G',
    'GGG' => 'G',
    'GGU' => 'G',
    'GUA' => 'V',
    'GUC' => 'V',
    'GUG' => 'V',
    'GUU' => 'V',
    'UAA' => '',
    'UAC' => 'Y',
    'UAG' => '',
    'UAU' => 'Y',
    'UCA' => 'S',
    'UCC' => 'S',
    'UCG' => 'S',
    'UCU' => 'S',
    'UGA' => '',
    'UGC' => 'C',
    'UGG' => 'W',
    'UGU' => 'C',
    'UUA' => 'L',
    'UUC' => 'F',
    'UUG' => 'L',
    'UUU' => 'F'
  }

  PROTEIN_CODON = {
    'A'=>['GCA', 'GCC', 'GCG', 'GCU'], 
    'C'=>['UGC', 'UGU'], 
    'D'=>['GAC', 'GAU'], 
    'E'=>['GAA', 'GAG'], 
    'F'=>['UUC', 'UUU'], 
    'G'=>['GGA', 'GGC', 'GGG', 'GGU'], 
    'H'=>['CAC', 'CAU'], 
    'I'=>['AUA', 'AUC', 'AUU'], 
    'K'=>['AAA', 'AAA', 'AAG'], 
    'L'=>['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'], 
    'M'=>['AUG'], 
    'N'=>['AAC', 'AAU'], 
    'P'=>['CCA', 'CCC', 'CCG', 'CCU'], 
    'Q'=>['CAA', 'CAG'], 
    'R'=>['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'], 
    'S'=>['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'], 
    'T'=>['ACA', 'ACC', 'ACG', 'ACU'], 
    'V'=>['GUA', 'GUC', 'GUG', 'GUU'], 
    'W'=>['UGG'], 
    'Y'=>['UAC', 'UAU']
  }
  
  MASS = {
    'G' => 57,
    'A' => 71,
    'S' => 87,
    'P' => 97,
    'V' => 99,
    'T' => 101,
    'C' => 103,
    'I' => 113,
    'L' => 113,
    'N' => 114,
    'D' => 115,
    'K' => 128,
    'Q' => 128,
    'E' => 129,
    'M' => 131,
    'H' => 137,
    'F' => 147,
    'R' => 156,
    'Y' => 163,
    'W' => 186,
  }

  MASS_ONLY = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]

  def self.rna_to_protein(rna)
    protein = []
    rna.split('').each_slice(3) {|kmer| protein << CODON_PROTEIN[kmer.join]}
    protein.join
  end

  def self.protein_to_rna(protein)
    product = PROTEIN_CODON[protein[0,1]]
    for i in 1..(protein.length-1)
      product = product.product( PROTEIN_CODON[protein[i,1]] )
    end
    product.map{|seq| seq.flatten.join}.uniq
  end

  def self.dna_to_rna(dna)
    dna.gsub('T','U')
  end

  def self.find_pattern_encodes_peptide(dna, peptide)
    res = []
    dna.split('').each_cons(peptide.length * 3) do |sub_dna|
      pattern = sub_dna.join
      pattern_rc = Chapter1.reverse_complement(pattern)

      res << pattern if rna_to_protein(dna_to_rna(pattern)) == peptide
      res << pattern if rna_to_protein(dna_to_rna(pattern_rc)) == peptide
    end
    res
  end

  def self.sub_peptides(peptide)
    subp = []
    for sub_len in 1..(peptide.length-1)
      for p in 0..(peptide.length-1)
        subp << (peptide+peptide)[p,sub_len]
      end
    end
    subp
  end

  def self.peptide_mass(peptide)
    mass = 0
    peptide.split('').each {|sp| mass += MASS[sp] }
    mass
  end

  def self.peptide_mass_i(peptide)
    peptide.inject(:+)
  end

  def self.cyclo_spectrum_dumb(peptide)
    masses = []
    sub_p = sub_peptides(peptide)
    masses << 0
    sub_p.each {|sp| masses << peptide_mass(sp)}
    masses << peptide_mass(peptide)
    masses.sort()
  end

  def self.linear_spectrum(peptide)
    prefix_mass = [0]
    for i in 0..(peptide.length-1)
      amino_acid = peptide[i,1]
      prefix_mass[i+1] = prefix_mass[i] + MASS[amino_acid]
    end
    linear_spectrum = [0]
    for i in 0..(peptide.length-1)
      for j in (i+1)..(peptide.length)
        linear_spectrum << prefix_mass[j] - prefix_mass[i]
      end
    end
    linear_spectrum.sort()
  end

  def self.cyclo_spectrum(peptide)
    prefix_mass = [0]
    for i in 0..(peptide.length-1)
      amino_acid_mass = peptide.is_a?(String) ? MASS[peptide[i,1]] : peptide[i]
      prefix_mass[i+1] = prefix_mass[i] + amino_acid_mass
    end
    peptide_mass = prefix_mass[peptide.length]
    cyclic_spectrum = [0]
    for i in 0..(peptide.length-1)
      for j in (i+1)..(peptide.length)
        diff_mass = prefix_mass[j] - prefix_mass[i]
        cyclic_spectrum << diff_mass
        cyclic_spectrum << peptide_mass - diff_mass if i > 0 && j < peptide.length
      end
    end
      cyclic_spectrum.sort()
  end

  def self.peptide_count(mass)
    for mi in 0..(MASS_ONLY.length-1)
      c_mass = 0

    end

  end

  def self.parent_mass(spectrum)
    spectrum.last
  end

  def self.cyclo_peptide_sequencing(spectrum)
    results = []
    peptides = (spectrum & Chapter2::MASS_ONLY).map{|aa| [aa]}

    peptides.each {|p| puts p.join('-')}

    while results.empty? && !peptides.empty? do

      #Branch
      peptides = expand(peptides, spectrum)

      # Bound
      peptides.delete_if do |peptide|
        delete = false

        if peptide_mass_i(peptide) == parent_mass(spectrum)
          cs = cyclo_spectrum(peptide)
          if cs == spectrum
            results << peptide
          end
          delete = true
        elsif !consistent(peptide, spectrum)
          delete = true
        end

        delete
      end

    end
    results
  end

  def self.expand(peptides, spectrum)
    ppts = []
    MASS_ONLY.each do |mass|

      if peptides.empty?
        ppts << [mass]
      else
        peptides.each do |origin_ppt|
          ppt = origin_ppt.dup
          ppt << mass
          ppts << ppt
        end
      end
    end
    ppts
  end

  def self.consistent(peptide, spectrum)
    mc_peptide  = map_count(peptide)
    mc_spectrum = map_count(spectrum)
    mc_peptide.each do |k,v|
      return false if mc_spectrum[k].nil? || mc_spectrum[k] < v
    end

    # subs = sub_peptides(peptide)
    # subs << peptide
    # subs.uniq.each do |s_p|
    #   unless spectrum.include?(peptide_mass_i(s_p))
    #     return false
    #   end
    # end

    return false unless spectrum.include?(peptide_mass_i(peptide))

    true
  end

  def self.map_count(arr)
    Hash[arr.sort.chunk{|i| i}.map{|k,v| [k,v.count]}]
  end

end