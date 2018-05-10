package org.panda.misc.analyses;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * @author Ozgun Babur
 */
public class MDACCTCGARPPA
{
	public static void main(String[] args) throws IOException
	{
		prepareFolders();
	}
	static void prepareFolders() throws IOException
	{
		String base = "/home/babur/Documents/RPPA/TCGA/correlation/";

		for (File file : new File(base).listFiles())
		{
			if (file.getName().endsWith(".csv"))
			{
				String study = file.getName().split("-")[1];

				String dir = base + study;
				new File(dir).mkdirs();

				int cols = Files.lines(Paths.get(file.getPath())).findFirst().get().split(",").length;
				boolean[] select = new boolean[cols];
				Arrays.fill(select, true);
				select[1] = false;
				select[2] = false;

				String dataFile = dir + "/data.txt";

				// write a transposed copy, also shorten sample names
				FileUtil.transpose(file.getPath(), ",", dataFile, "\t", t ->
				{
					if (t[0].endsWith("_ID")) t[0] = "ID";
					else if (t[0].length() > 12) t[0] = t[0].substring(0, 12);
				}, select);

				BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "/parameters.txt"));
				writer.write(PARAMS);

				String[] c = Files.lines(Paths.get(dataFile)).findFirst().get().split("\t");
				for (int i = 1; i < c.length; i++)
				{
					writer.write("\nvalue-column = " + c[i]);
				}

				writer.close();
			}
		}
	}

	public static final String PARAMS = "proteomics-platform-file = ../tcga-platform.txt\n" +
		"proteomics-values-file = data.txt\n" +
		"id-column = ID\n" +
		"symbols-column = Symbols\n" +
		"sites-column = Sites\n" +
		"effect-column = Effect\n" +
		"\n" +
		"value-transformation = correlation\n" +
		"fdr-threshold-for-correlation = 0.1\n" +
		"\n" +
		"calculate-network-significance = false\n" +
		"generate-data-centric-graph = true\n";
}
