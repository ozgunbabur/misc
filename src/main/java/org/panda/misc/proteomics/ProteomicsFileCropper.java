package org.panda.misc.proteomics;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ProteomicsFileCropper
{
	public static void crop(String inFile, String outFile, Set<String> cols) throws IOException
	{
		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write(header[0]);

		for (int i = 1; i < header.length; i++)
		{
			if (i < 4 || cols.contains(header[i]))
			{
				writer.write("\t" + header[i]);
			}
		}

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			FileUtil.lnwrite(t[0], writer);

			for (int i = 1; i < header.length; i++)
			{
				if (i < 4 || cols.contains(header[i]))
				{
					FileUtil.tab_write(t[i], writer);
				}
			}
		});

		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/";
		Set<String> samples = Files.lines(Paths.get(dir + "samples-77.txt")).collect(Collectors.toSet());
		crop(dir + "CPTAC-TCGA-BRCA-data.txt", dir + "CPTAC-TCGA-BRCA-data-77.txt", samples);
	}
}
