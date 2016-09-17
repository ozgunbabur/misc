package org.panda.misc.pancan;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Extracts part of the MAF file to generate a lollipop plot using
 * http://www.cbioportal.org/mutation_mapper.jsp
 *
 * @author Ozgun Babur
 */
public class LollipopFromMAF
{
	String maffile;
	String outFile;
	Set<String> genes;

	public LollipopFromMAF(String maffile, String outFile, Set<String> genes)
	{
		this.maffile = maffile;
		this.outFile = outFile;
		this.genes = genes;
	}

	public void run() throws IOException
	{
		String[] header = Files.lines(Paths.get(maffile)).filter(l -> !l.startsWith("#")).findFirst().get().split("\t");

		int geneInd = ArrayUtil.indexOf(header, "Hugo_Symbol");
		int changeInd = ArrayUtil.indexOf(header, "Protein_Change", "HGVSp_Short");
		int sampleInd = ArrayUtil.indexOf(header, "Sample_ID", "Tumor_Sample_Barcode");
		int typeInd = ArrayUtil.indexOf(header, "Mutation_Type", "Variant_Classification");

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("Hugo_Symbol\tProtein_Change\tSample_ID\tMutation_Type");
		Files.lines(Paths.get(maffile)).map(l -> l.split("\t")).filter(t -> genes.contains(t[geneInd])).forEach(t ->
			FileUtil.lnwrite(t[geneInd] + "\t" + t[changeInd] + "\t" + t[sampleInd] + "\t" + t[typeInd], writer));
		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		LollipopFromMAF lol = new LollipopFromMAF("/home/babur/Documents/TCGA/PanCan/mutation.maf",
			"/home/babur/Documents/PanCan/lollipop/MTOR-2.maf", new HashSet<>(Arrays.asList("MTOR")));

		lol.run();
	}
}
