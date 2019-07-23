package org.panda.misc.altmatrix;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;

/**
 * Gets sub-matrices of the given alteration matrix.
 *
 * @author Ozgun Babur
 */
public class AltMatrixCropper
{
	public static void crop(String matrixFile, Set<String> genes, Set<String> samples, String outFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(matrixFile)).findFirst().get().split("\t");

		boolean[] select = new boolean[header.length];

		for (int i = 0; i < select.length; i++)
		{
			select[i] = samples == null || samples.contains(header[i]);
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("Genes");

		for (int i = 1; i < header.length; i++)
		{
			if (select[i]) writer.write("\t" + header[i]);
		}

		Files.lines(Paths.get(matrixFile)).skip(1).map(l -> l.split("\t"))
			.filter(t -> genes == null || genes.contains(t[0])).forEach(t ->
		{
			FileUtil.lnwrite(t[0], writer);
			for (int i = 1; i < t.length; i++)
			{
				if (select[i]) FileUtil.tab_write(t[i], writer);
			}
		});

		writer.close();
	}



	public static void main(String[] args) throws IOException
	{
		crop("/media/babur/6TB1/TCGA-pancan/whole/TP53-neighborhood/DataMatrix-TP53-neighborhood-with-extra-samples.txt", null,
			FileUtil.getTermsInTabDelimitedColumn("/media/babur/6TB1/TCGA-pancan/sample-to-tissue-mapping.txt", 0, 1),
			"/media/babur/6TB1/TCGA-pancan/whole/TP53-neighborhood/DataMatrix-TP53-neighborhood.txt");
	}
}
