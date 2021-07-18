package org.panda.misc.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Damien
{
	public static final String BASE = "/Users/ozgun/Documents/Analyses/Damien-UofO/";

	public static void main(String[] args) throws IOException
	{
		convertPhospho();
	}

	public static void convertPhospho() throws IOException
	{
		String filename = BASE + "data.tsv";
		String[] header = FileUtil.lines(filename).skip(4).findFirst().get().split("\t");
		int fcInd = ArrayUtil.indexOf(header, "logFC");
		int pInd = ArrayUtil.indexOf(header, "FDR");
		int symInd = ArrayUtil.indexOf(header, "Protein Description");
		int siteInd = ArrayUtil.indexOf(header, "Modification Location (Protein)");

		BufferedWriter writer = FileUtil.newBufferedWriter(BASE + "data.txt");
		writer.write("ID\tSymbols\tSites\tEffect\tSignedP");

		Set<String> mem = new HashSet<>();
 		FileUtil.lines(filename).skip(5).map(l -> l.split("\t")).forEach(t ->
		{
			int gStart = t[symInd].indexOf(" GN=");
			if (gStart < 0) return;
			String sym = t[symInd].substring(gStart + 4, t[symInd].indexOf(" ", gStart + 4));

			List<String> sites = new ArrayList<>();
			for (String s : t[siteInd].split("; "))
			{
				if (s.endsWith("(Phospho)"))
				{
					sites.add(s.substring(0, s.indexOf("(")));
				}
			}

			String id = sym + "-" + CollectionUtil.merge(sites, "-");
			String idd = id;
			int rep = 1;
			while (mem.contains(idd)) idd = id + "-rep" + (rep++);
			id = idd;
			mem.add(id);

			double p = Double.valueOf(t[pInd]);
			if (t[fcInd].startsWith("-")) p = -p;

			FileUtil.lnwrite(id + '\t' + sym + "\t" + CollectionUtil.merge(sites, "|") + "\t\t" + p, writer);
		});

 		writer.close();
	}
}
