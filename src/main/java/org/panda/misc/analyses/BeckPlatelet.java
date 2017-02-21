package org.panda.misc.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * For the analysis of the data int he paper:
 * https://www.ncbi.nlm.nih.gov/pubmed/28060719
 *
 * @author Ozgun Babur
 */
public class BeckPlatelet
{
	public static final String DIR = "/home/babur/Documents/Analyses/BeckProteomics/";

	public static void main(String[] args) throws IOException
	{
		prepareDataFile();
	}

	public static void prepareDataFile() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "data.txt"));
		writer.write("ID\tSymbols\tSites\tEffect\t10s\t20s\t30s");
		Files.lines(Paths.get(DIR + "TableS3.txt")).skip(2).map(l -> l.split("\t")).filter(t -> !t[2].startsWith("n"))
			.forEach(t ->
		{
			String sym = t[2];
			String[] sites = t[4].split(" or | and | and/or ");
			Double[] v = new Double[]{Double.valueOf(t[5]), Double.valueOf(t[8]), Double.valueOf(t[11])};
			for (int i = 0; i < v.length; i++)
			{
				if (v[i] < 1) v[i] = -(1 / v[i]);
			}

			List<String> siteList = new ArrayList<>();
			for (String site : sites)
			{
				siteList.add("S" + site);
				siteList.add("T" + site);
				siteList.add("Y" + site);
			}

			String id = sym;
			for (String site : sites)
			{
				id += "_" + site;
			}

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(siteList, "|") + "\t\t" +
				ArrayUtil.getString("\t", v), writer);

		});

		writer.close();
	}
}
