package org.panda.misc.analyses;

import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.resource.UniProtSequence;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PaulBPhosphoproteomics
{
	public static final String DIR = "/home/babur/Documents/Analyses/PaulBarnes/";

	public static void main(String[] args) throws IOException
	{
		convertMouseToHuman();
	}

	public static void convertMouseToHuman() throws IOException
	{
		int tnum = 3;
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "T" + tnum + "/data.txt"));
		writer.write("ID\tGenes\tSites\tEffect\tFold-change");

		Files.lines(Paths.get(DIR + "Treatment" + tnum + ".csv")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String pep = t[0];
			String markedPep = t[1];
			String name = t[2];
			double fc = Double.valueOf(t[3]);

			int start = getStartPosition(pep, name);
			if (start < 1) System.err.println("Cannot find peptide: " + name);

			List<String> sites = getModifications(markedPep, start);
			String upID = UniProtSequence.get().getIDOfName(name);

			List<String> hSites = new ArrayList<>();
			for (String site : sites)
			{
				String hSite = SiteMappingMouseToHuman.get().map(upID, site);
				if (hSite != null) hSites.add(hSite);
			}
			if (!hSites.isEmpty())
			{
				String symbol = hSites.iterator().next().split("-")[0];
				String siteStr = collectSites(hSites);

				FileUtil.lnwrite(symbol + "-" + siteStr.replaceAll("\\|", "-") + "\t" + symbol + "\t" + siteStr + "\t\t" + fc,
					writer);
			}
		});

		writer.close();
	}

	private static String collectSites(Collection<String> genesWithSites)
	{
		List<String> sites = genesWithSites.stream().map(s -> s.split("-")[1]).collect(Collectors.toList());
		return CollectionUtil.merge(sites, "|");
	}

	static int getStartPosition(String peptide, String uniprotName)
	{
		return UniProtSequence.get().getStartLocation(uniprotName, peptide);
	}

	static List<String> getModifications(String markedPeptide, int startPos)
	{
		markedPeptide = markedPeptide.substring(2, markedPeptide.length() - 2);
		List<String> sites = new ArrayList<>();
		String[] t = markedPeptide.split("\\*");
		for (int i = 0; i < t.length - 1; i++)
		{
			char lett = t[i].charAt(t[i].length() - 1);
			int pos = startPos + t[i].length() - 1;

			sites.add(lett + "" + pos);

			startPos += t[i].length();
		}

		return sites;
	}


}
