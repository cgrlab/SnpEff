package ca.mcgill.mcb.pcingola.genotypes;

import java.io.Serializable;
import java.util.List;

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

	int size; // Size in elements (genotypes)
	byte genotype[];
	boolean haploid; // True if haploid, false if diploid

	public GenotypeVector(int size) {
		init(size);
	}

	public GenotypeVector(VcfEntry ve) {
		set(ve);
	}

	public String get(int sampleNum) {
		byte code = genotype[sampleNum];

		if (haploid) {
			// Haploid chromosome
			if ((code & MISSING_MASK) != 0) return ".";
			return (code & GT_MASK) == 0 ? "0" : "1";
		}

		// Diploid
		String slash = ((code & PHASED_MASK) != 0) ? "|" : "/";
		if ((code & MISSING_MASK) != 0) return "." + slash + ".";
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

	public int getCode(int sampleNum) {
		return genotype[sampleNum];
	}

	protected void init(int size) {
		this.size = size;
		genotype = new byte[size];
	}

	public boolean isMissing(int sampleNum) {
		return (genotype[sampleNum] & MISSING_MASK) != 0;
	}

	public boolean isPhased(int sampleNum) {
		return (genotype[sampleNum] & PHASED_MASK) != 0;
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
	public void set(int sampleNum, VcfGenotype vg) {
		byte code = (byte) vg.getGenotypeCode();
		if (code < 0) code = MISSING_MASK;
		if (vg.isPhased()) code |= PHASED_MASK;
		set(sampleNum, code);
	}

	public void set(int sampleNum, VcfGenotype vg, String alt) {
		byte code = (byte) vg.getGenotypeCode(alt);
		if (code < 0) code = MISSING_MASK;
		if (vg.isPhased()) code |= PHASED_MASK;
		set(sampleNum, code);
	}

	public void set(VcfEntry ve) {
		List<VcfGenotype> gts = ve.getVcfGenotypes();

		init(gts.size());

		int i = 0;
		int ploidy = 0;
		for (VcfGenotype gt : gts) {
			set(i++, gt);
			ploidy = Math.max(ploidy, gt.plodity());
		}

		setPloidy(ploidy);
	}

	public void set(VcfEntry ve, String alt) {
		List<VcfGenotype> gts = ve.getVcfGenotypes();

		init(gts.size());

		int i = 0;
		for (VcfGenotype gt : gts)
			set(i++, gt, alt);
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
		return size;
	}
}
