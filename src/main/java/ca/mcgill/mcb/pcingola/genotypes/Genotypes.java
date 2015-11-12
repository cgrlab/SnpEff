package ca.mcgill.mcb.pcingola.genotypes;

import java.io.Serializable;
import java.util.Collection;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.Haplotype;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;

/**
 * Genotypes for a set of (diploid/haploid) samples
 *
 * Each genotype is stored in 1 byte:
 *   Phasing information is stored in bit 7.
 *   Missing information is stored in bit 6.
 *   Genotypes 0/0, 0/1, 1/0, 1/1 are stored bits 0 and 1
 *
 *   Coding: 76543210
 *           ^^^   ^^
 *           ||   Genotype: 0=0/0, 1=0/1, 2=1/0, 3=1/1
 *           ||
 *           |Missing: 1 if missing, 0 otherwise
 *           |
 *           Phase: 1 if phased, 0 otherwise
 *
 * @author pcingola
 */
public class Genotypes implements Serializable {

	protected static final byte PHASED_MASK = (byte) 0x80;
	protected static final byte MISSING_MASK = (byte) 0x40;
	protected static final byte GT_MASK = (byte) 0x03;

	protected static final int HAPLOID = 1;
	protected static final int DIPLOID = 2;

	private static final long serialVersionUID = 4734574592894281057L;

	byte genotype[];
	boolean haploid; // True if haploid, false if diploid
	Haplotype haplotype;

	public Genotypes(int size) {
		init(size);
	}

	public Genotypes(VcfEntry ve) {
		set(ve);
	}

	public Genotypes(VcfEntry ve, String alt) {
		set(ve, alt);
	}

	/**
	 * Are all genotypes 'REF'?
	 */
	public boolean allRef() {
		int size = size();
		for (int i = 0; i < size; i++)
			if (!isHomozygousRef(genotype[i])) return false;

		return true;
	}

	/**
	 * Calculate the resulting genotype when combining with 
	 * the previous entry 'gtOld'.
	 */
	protected Genotypes calcGenotypeVector(Genotypes gtOld) {
		if (gtOld.getHaplotype().getLastStart() > getHaplotype().getFirstStart()) //
			throw new RuntimeException("Haplotype sort order reversed!" //
					+ "\n\tthis  :" + getHaplotype() //
					+ "\n\tgtOld :" + gtOld.getHaplotype() //
		);

		// Create new genotypes
		int size = size();
		Genotypes newGs = new Genotypes(size);

		// Update haplotype
		Haplotype newHap = newGs.getHaplotype();
		newHap.add(this.getHaplotype());
		newHap.add(gtOld.getHaplotype());

		if (!isHaploid() && !gtOld.isHaploid()) {
			newGs.setPloidy(DIPLOID);

			// Both are diploid genomes
			for (int i = 0; i < size; i++) {
				byte gtCode0 = gtOld.getCode(i);
				byte gtCode1 = getCode(i);

				if (isPhased(gtCode0) && isPhased(gtCode1)) {
					// At least one bit remains the same in both genotypes
					// (i.e. one genotype is the same)
					newGs.set(i, (byte) (PHASED_MASK | (gtCode0 & gtCode1)));
				} else if (!isPhased(gtCode0) && isPhased(gtCode1)) {
					// Old genotype is not phased, but this new one is phased
					// This is a case of a new 'phasing block' starting (see GATK's 
					// Read Backed Phasing for details). So this results in a 
					// phased genotype 
					// Reference: http://gatkforums.broadinstitute.org/discussion/45/purpose-and-operation-of-read-backed-phasing
					newGs.set(i, (byte) (PHASED_MASK | (gtCode0 & gtCode1)));
				} else {
					// We need the same condition as in 'phased', but also
					// at least one of them to be homozygous-ALT
					if (isHomozygousAlt(gtCode0) || isHomozygousAlt(gtCode1)) newGs.set(i, (byte) (gtCode0 & gtCode1));
				}
			}
		} else if (isHaploid() && gtOld.isHaploid()) {
			newGs.setPloidy(HAPLOID);

			// Both are haploid genomes: Any matching non-ref is implicitly
			// phased since these are haploid chromosomes
			for (int i = 0; i < size; i++) {
				int gt0 = getCode(i) & GT_MASK;
				int gt1 = gtOld.getCode(i) & GT_MASK;
				if (gt0 != 0 && gt1 != 0) newGs.set(i, (byte) 1);
			}
		} else {
			// One is haploid and the other is diploid
			// This is a weird case.
			// Also, may be a this implies a loss of heterozygosity
			// I'm not sure if this is compliant with the VCF spec.)
			throw new RuntimeException("Genotype verctors have different zygosity?" //
					+ "\n\t" + this //
					+ "\n\t" + gtOld //
			);
		}

		return newGs;
	}

	/**
	 * Compare all genotypes
	 */
	public boolean equalsGenotypes(Genotypes gt) {
		int size = size();
		for (int i = 0; i < size; i++)
			if ((genotype[i] & GT_MASK) != (gt.genotype[i] & GT_MASK)) {
				Gpr.debug("NOT EQUALS\n\tthis: " + this + "\n\tgt  :" + gt);
				return false;
			}

		Gpr.debug("EQUALS\n\tthis: " + this + "\n\tgt  :" + gt);
		return true;
	}

	/**
	 * Find index of 'alt' entry in VcfEntry
	 */
	int findAltIndex(VcfGenotype gt, String alt) {
		VcfEntry vcfEntry = gt.getVcfEntry();
		String alts[] = vcfEntry.getAlts();
		int altIndex = -1;
		for (int i = 0; i < alts.length; i++)
			if (alts[i].equalsIgnoreCase(alt)) {
				altIndex = i;
				break;
			}

		// Sanity check
		if (altIndex < 0) throw new RuntimeException("ALT '" + alt + "' does not match any ALT in VCF entry:\t" + vcfEntry);
		return altIndex;
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
			return "1" + slash + "0";

		case 2:
			return "0" + slash + "1";

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
		int altCode = findAltIndex(gt, alt) + 1; // Code 0 means 'REF', so we have to add 1

		// Calculate code: Genotype matches 'alt'? Then set bit
		// Note: We only support haploid / diploid entries
		byte code = 0;
		if (gtCodes[0] == altCode) code |= 1;
		if (gtCodes.length > 1 && gtCodes[1] == altCode) code |= 2;

		return code;
	}

	public Haplotype getHaplotype() {
		return haplotype;
	}

	/**
	 * Can this genotype conform a haplotype with all the genotypes in 'gts'?
	 * @return The resulting genotype on success, null otherwise
	 */
	public Genotypes haplotype(Collection<Genotypes> gts) {
		Genotypes gtResult = this;

		for (Genotypes gt : gts) {
			gtResult = gtResult.calcGenotypeVector(gt);
			Gpr.debug("Calculated genotype: " + gtResult);
			if (gtResult.allRef()) return null; // All entries are 'REF'? There is nothing else to do
		}

		return gtResult;
	}

	/**
	 * Can this genotype conform a haplotype with genotype 'gt' ?
	 * @return The resulting genotype on success, null otherwise
	 */
	public Genotypes haplotype(Genotypes gt) {
		Genotypes gtResult = this.calcGenotypeVector(gt);
		if (gtResult.allRef()) return null; // All entries are 'REF'? There is nothing else to do
		return gtResult;
	}

	protected void init(int size) {
		genotype = new byte[size];
		haplotype = new Haplotype();
	}

	public boolean isHaploid() {
		return haploid;
	}

	boolean isHomozygousAlt(byte gtCode) {
		return (gtCode & GT_MASK) == GT_MASK;
	}

	boolean isHomozygousRef(byte gtCode) {
		return (gtCode & GT_MASK) == 0;
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

		// Initialize
		List<VcfGenotype> gts = ve.getVcfGenotypes();
		init(gts.size());
		haplotype.add(ve.variants().get(0));

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
		Variant variant = null;
		for (Variant var : ve.variants())
			if (var.getGenotype().equalsIgnoreCase(alt)) {
				variant = var;
				break;
			}
		if (variant == null) throw new RuntimeException("Cannot find corresponding variant for ALT='" + alt + "'");

		// Initialize
		List<VcfGenotype> gts = ve.getVcfGenotypes();
		init(gts.size());
		haplotype.add(variant);

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
		if (haplotype == null) {
			sb.append("null");
		} else if (haplotype.size() == 1) {
			sb.append("Variant: ");
			Variant variant = haplotype.getVariants().get(0);
			sb.append(variant.getChromosomeName() + ":" + (variant.getStart() + 1) + "_" + variant.getReference() + ">" + variant.getAlt());
		} else {
			sb.append("Haplotype: " + haplotype);
		}

		sb.append(", " + (haploid ? "haploid" : "diploid"));
		sb.append(", size:" + size);
		sb.append(", genotypes:");

		for (int i = 0; i < size; i++)
			sb.append(" " + get(i));

		return sb.toString();
	}
}
