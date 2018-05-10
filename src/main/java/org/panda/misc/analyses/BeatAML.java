package org.panda.misc.analyses;

import org.panda.utility.FileUtil;

/**
 * @author Ozgun Babur
 */
public class BeatAML
{
	public static final String DIR = "/home/babur/Documents/Analyses/BeatAML/";


	public static void main(String[] args)
	{
		prepareMatrices();
	}

	public static void prepareMatrices()
	{
		String outFile = DIR + "mutex-all/matrix.txt";
		FileUtil.fixRStyleHeader(DIR + "curated_variants_idxhugo_all.tsv", "\t", outFile);
		FileUtil.transpose(outFile, "\t", outFile, "\t", null, null);

		outFile = DIR + "mutex-xy-filtered/matrix.txt";
		FileUtil.fixRStyleHeader(DIR + "curated_variants_idxhugo_filt.tsv", "\t", outFile);
		FileUtil.transpose(outFile, "\t", outFile, "\t", null, null);
	}
}
