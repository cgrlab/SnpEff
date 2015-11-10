package ca.mcgill.mcb.pcingola.genotypes;

import java.io.Serializable;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;

/**
 * A vector of genotypes in a 'compact' structure
 *
 * Phasing information is stored in bit 7.
 * Missing information is stored in bit 6.
 * Genotypes 0/0, 0/1, 1/0, 1/1 are stored bits 0 and 1
 *
 * Coding: 76543210
 *         ^^^   ^^
 *         ||   Genotype: 0=0/0, 1=0/1, 2=1/0, 3=1/1
 *         ||
 *         |Missing: 1 if missing, 0 otherwise
 *         |
 *         Phase: 1 if phased, 0 otherwise
 *
 * @author pcingola
 */
public class GenotypeVector implements Serializable {

	protected static final byte PHASED_MASK = (byte) 0x80;
	protected static final byte MISSING_MASK = (byte) 0x40;
	protected static final byte GT_MASK = (byte) 0x03;

	private static final long serialVersionUID = 4734574592894281057L;

	byte genotype[];
	boolean haploid; // True if haploid, false if diploid
	Variant variant;

	public GenotypeVector(int size) {
		init(size);
	}

	public GenotypeVector(VcfEntry ve) {
		set(ve);
	}

	public GenotypeVector(VcfEntry ve, String alt) {
		set(ve, alt);
	}

	/**
	 * Find index of 'alt' entry in VcfEntry
	 */
	int findAltIndex(VcfGenotype gt, String alt) {
		VcfEntry vcfEntry = gt.getVcfEntry();
		String alts[] = vcfEntry.getAlts();
		int altCode = -1;
		for (int i = 0; i < alts.length; i++)
			if (alts[i].equalsIgnoreCase(alt)) {
				altCode = i;
				break;
			}

		// Sanity check
		if (altCode < 0) throw new RuntimeException("ALT '" + alt + "' does not match any ALT in VCF entry:\t" + vcfEntry);
		return altCode;
	}

	public String get(int sampleNum) {
		byte code = genotype[sampleNum];

		if (haploid) {
			// Haploid chromosome
			if (isMissing(code)) return ".";
			return (code & GT_MASK) == 0 ? "0" : "1";
		}

		// Diploid
		String slash = isPhased(code) ? "|" : "/";
		if (isMissing(code)) return "." + slash + ".";
		switch (code & GT_MASK) {
		case 0:
			return "0" + slash + "0";

		case 1:
			return "0" + slash + "1";

		case 2:
			return "1" + slash + "0";

		case 3:
			return "1" + slash + "1";

		default:
			throw new RuntimeException("Unknown code '" + code + "'");
		}
	}

	public byte getCode(int sampleNum) {
		return genotype[sampleNum];
	}

	/**
	 * Calculate genotype code.
	 */
	byte getGenotypeCode(VcfGenotype gt) {
		if (gt.isMissing()) return -1;

		int gtCodes[] = gt.getGenotype();
		if (gtCodes == null) return -1;

		// Calculate code: Genotype matches is not 'REF'? Then set bit
		// Note: We only support haploid / diploid entries
		byte code = 0;
		if (gtCodes[0] != 0) code |= 1;
		if (gtCodes.length > 1 && gtCodes[1] != 0) code |= 2;

		return code;
	}

	/**
	 * Calculate genotype code.
	 * Only assume "ALT" if it matches the provided 'alt'
	 */
	byte getGenotypeCode(VcfGenotype gt, String alt) {
		if (gt.isMissing()) return -1;

		int gtCodes[] = gt.getGenotype();
		if (gtCodes == null) return -1;

		// Find corresponding ALT code
		int altIndex = findAltIndex(gt, alt);

		// Calculate code: Genotype matches 'alt'? Then set bit
		// Note: We only support haploid / diploid entries
		byte code = 0;
		if (gtCodes[0] == altIndex) code |= 1;
		if (gtCodes.length > 1 && gtCodes[1] == altIndex) code |= 2;

		return code;
	}

	public Variant getVariant() {
		return variant;
	}

	/**
	 * Is there a matching haplotype with 'gv'
	 */
	public boolean hasHaplotype(GenotypeVector gv) {
		int size = size();

		if (!isHaploid() && !gv.isHaploid()) {
			// Both are diploid genomes
			for (int i = 0; i < size; i++) {
				byte gt0 = getCode(i);
				byte gt1 = gv.getCode(i);

				if (isPhased(gt0) && isPhased(gt1)) {
					// At least one bit remains the same in both genotypes
					// (i.e. one genotype is the same)
					// Note: We implicitly use the fact that missing
					// genotypes are set to zero.
					if ((gt0 & gt1 & GT_MASK) != 0) return true;
				} else {
					// We need the same condition as in 'phased', but also
					// at least one of them to be homozygous-ALT
					if ((gt0 & gt1 & GT_MASK) != 0 //
							&& (isHomozygousAlt(gt0) || isHomozygousAlt(gt1)) //
					) return true;
				}
			}
		} else if (isHaploid() && gv.isHaploid()) {
			// Both are haploid genomes: Any matching non-ref is implicitly
			// phased since these are haploid chromosomes
			for (int i = 0; i < size; i++) {
				int gt0 = getCode(i) & GT_MASK;
				int gt1 = gv.getCode(i) & GT_MASK;
				if (gt0 != 0 && gt1 != 0) return true;
			}

		} else {
			// One is haploid and the other is diploid
			// This is a weird case.
			// Also, may be a this implies a loss of heterozygosity
			// I'm not sure if this is compliant with the VCF spec.)
			throw new RuntimeException("Genotype verctors have different zygosity?" //
					+ "\n\t" + this //
					+ "\n\t" + gv //
			);
		}

		return false;
	}

	protected void init(int size) {
		genotype = new byte[size];
	}

	public boolean isHaploid() {
		return haploid;
	}

	boolean isHomozygousAlt(byte gtCode) {
		return (gtCode & GT_MASK) == GT_MASK;
	}

	boolean isMissing(byte gtCode) {
		return (gtCode & MISSING_MASK) != 0;
	}

	boolean isPhased(byte gtCode) {
		return (gtCode & PHASED_MASK) != 0;
	}

	/**
	 * Set genotype code
	 * Codes {0, 1, 2, 3} => Genotypes { 0/0, 0/1, 1/0, 1/1 }
	 */
	public void set(int sampleNum, byte code) {
		genotype[sampleNum] = code;
	}

	/**
	 * Set genotype
	 */
	void set(int sampleNum, VcfGenotype vg) {
		byte code = getGenotypeCode(vg);
		if (code < 0) code = MISSING_MASK; // Note that we set the genotype part to zero.
		if (vg.isPhased()) code |= PHASED_MASK;
		set(sampleNum, code);
	}

	void set(int sampleNum, VcfGenotype vg, String alt) {
		byte code = getGenotypeCode(vg, alt);
		if (code < 0) code = MISSING_MASK; // Note that we set the genotype part to zero.
		if (vg.isPhased()) code |= PHASED_MASK;
		set(sampleNum, code);
	}

	void set(VcfEntry ve) {
		if (ve.isMultiallelic()) throw new RuntimeException("Cannot add ulti-allelic VCF entries without specifiying an 'ALT'");

		// Set variant
		variant = ve.variants().get(0);

		// Initialize
		List<VcfGenotype> gts = ve.getVcfGenotypes();
		init(gts.size());

		// Set genotypes
		int i = 0;
		int ploidy = 0;
		for (VcfGenotype gt : gts) {
			set(i++, gt);
			ploidy = Math.max(ploidy, gt.plodity());
		}

		setPloidy(ploidy);
	}

	public void set(VcfEntry ve, String alt) {
		// Find corresponding variant
		for (Variant var : ve.variants())
			if (var.getGenotype().equalsIgnoreCase(alt)) {
				variant = var;
				break;
			}
		if (variant == null) throw new RuntimeException("Cannot find corresponding variant for ALT='" + alt + "'");

		// Initialize
		List<VcfGenotype> gts = ve.getVcfGenotypes();
		init(gts.size());

		// Set genotypes
		int i = 0;
		int ploidy = 0;
		for (VcfGenotype gt : gts) {
			set(i++, gt, alt);
			ploidy = Math.max(ploidy, gt.plodity());
		}

		setPloidy(ploidy);
	}

	protected void setPloidy(int ploidy) {
		switch (ploidy) {
		case 1:
			haploid = true;
			break;

		case 0: // Not found, assume diploid
		case 2:
			haploid = false;
			break;

		default:
			throw new RuntimeException("Polyploidy not supported: ploidy=" + ploidy);
		}
	}

	public int size() {
		return genotype == null ? 0 : genotype.length;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		int size = size();
		sb.append("Variant: ");
		if (variant != null) sb.append(variant.getChromosomeName() + ":" + (variant.getStart() + 1) + "_" + variant.getReference() + ">" + variant.getAlt());
		else sb.append("null");

		sb.append(", " + (haploid ? "haploid" : "diploid"));
		sb.append(", size:" + size);
		sb.append(", genotypes:");

		for (int i = 0; i < size; i++)
			sb.append(" " + get(i));

		return sb.toString();
	}
}
