package org.panda.misc.analyses;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class Mundt
{
	public static final String DIR = "/home/babur/Documents/Analyses/Mundt/";
	public static final double FDR_THR = 0.1;

	public static void main(String[] args) throws IOException
	{
		convert2H();
	}

	static void convert2H() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "50H/data-fdr0.1.txt"));
		writer.write("ID\tGenes\tSites\tEffect");

		String filename = "Phosphoproteome-50h.csv";

		String[] header = Files.lines(Paths.get(DIR + filename)).findFirst().get().split("\t");

		for (int i = 3; i < header.length; i++)
		{
			writer.write("\t" + header[i]);
		}

		Files.lines(Paths.get(DIR + filename)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (Double.valueOf(t[2]) > FDR_THR) return;

			String gene = t[0];
			String sites = CollectionUtil.merge(getSites(t[1]), "|");
			FileUtil.lnwrite(gene + "-" + sites.replaceAll("\\|", "-") + "\t" + gene + "\t" + sites + "\t", writer);

			for (int i = 3; i < t.length; i++)
			{
				FileUtil.write("\t" + t[i].replaceAll("NA", "NaN"), writer);
			}
		});

		filename = "Proteome-50h.csv";

		Files.lines(Paths.get(DIR + filename)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (Double.valueOf(t[1]) > FDR_THR) return;

			String gene = t[0];
			FileUtil.lnwrite(gene + "\t" + gene + "\t\t", writer);

			for (int i = 2; i < t.length; i++)
			{
				FileUtil.write("\t" + t[i].replaceAll("NA", "NaN"), writer);
			}
		});

		writer.close();
	}

	static List<String> getSites(String s)
	{
		s = s.split("_")[2];
		String[] t = s.split("s|t|y");

		List<String> sites = new ArrayList<>();
		for (int i = 0; i < t.length; i++)
		{
			if (t[i].length() > 1) sites.add(t[i]);
		}
		return sites;
	}
}
