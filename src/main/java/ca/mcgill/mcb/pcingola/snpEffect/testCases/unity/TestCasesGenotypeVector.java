package ca.mcgill.mcb.pcingola.snpEffect.testCases.unity;

import java.util.Random;

import org.junit.Test;

import ca.mcgill.mcb.pcingola.genotypes.GenotypeVector;
import ca.mcgill.mcb.pcingola.util.Gpr;
import junit.framework.Assert;

/**
 * Test cases for GenotypeVector class
 *
 * @author pcingola
 */
public class TestCasesGenotypeVector {

	boolean verbose = false;

	@Test
	public void test_01() {
		Gpr.debug("Test");
		int len = 100;
		for (int code = 0; code < 4; code++) {
			GenotypeVector gv = new GenotypeVector(len);

			for (int i = 0; i < len; i++)
				gv.set(i, (byte) code);

			for (int i = 0; i < len; i++)
				Assert.assertEquals(code, gv.getCode(i));
		}
	}

	@Test
	public void test_02() {
		Gpr.debug("Test");
		Random rand = new Random(20121221);
		GenotypeVector gv = new GenotypeVector(1000);

		// Create random codes
		int codes[] = new int[gv.size()];
		for (int i = 0; i < gv.size(); i++) {
			byte code = (byte) rand.nextInt(4);
			codes[i] = code;
			gv.set(i, code);
		}

		// Check that codes are stored OK
		for (int i = 0; i < gv.size(); i++) {
			Assert.assertEquals(codes[i], gv.getCode(i));
		}
	}
}
