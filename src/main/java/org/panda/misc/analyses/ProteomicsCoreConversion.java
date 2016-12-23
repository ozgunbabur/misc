package org.panda.misc.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class ProteomicsCoreConversion
{

	static void convert(String inFile, String outFile, double qvalThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(inFile)).skip(4).findFirst().get().split("\t");

		int siteInd = ArrayUtil.indexOf(header, "Mod Protein Location");
		if (siteInd < 0) ArrayUtil.indexOf(header, "Protein Modification Location");
		int qvalInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "fold change");
		int logfcInd = ArrayUtil.indexOf(header, "logFC");
		int descInd = ArrayUtil.indexOf(header, "Protein Descriptions");

		Map<String, Holder> map = new HashMap<>();

		Files.lines(Paths.get(inFile)).skip(5).map(l -> l.split("\t"))
			.filter(t -> Double.parseDouble(t[qvalInd]) < qvalThr).forEach(t ->
		{
			String sym = t[descInd];

			int start = sym.indexOf("GN=") + 3;
			sym = sym.substring(start, sym.indexOf(" ", start));

			String[] sites = t[siteInd].split("; ");
			double fc = Double.parseDouble(t[fcInd]);
			if (t[logfcInd].startsWith("-")) fc = -fc;

			String id = sym;
			for (String site : sites) id += "_" + site;

			Holder h = new Holder(id, sym, fc);

			for (String site : sites)
			{
				String key = sym + " " + site;
				if (!map.containsKey(key) || Math.abs(map.get(key).fc) < Math.abs(fc))
				{
					map.put(key, h);
				}
			}
		});

		map.keySet().forEach(k -> map.get(k).sites.add(k.substring(k.indexOf(" ") + 1)));

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tFold change");

		map.values().stream().distinct().forEach(h -> FileUtil.lnwrite(h.toString(), writer));

		writer.close();
	}

	static class Holder
	{
		String sym;
		String id;
		double fc;
		List<String> sites;

		Holder(String id, String sym, double fc)
		{
			this.sym = sym;
			this.id = id;
			this.fc = fc;
			sites = new ArrayList<>();
		}

		@Override
		public String toString()
		{
			return id + "\t" + sym + "\t" + ArrayUtil.getString("|", sites.toArray(new String[sites.size()])) +
				"\t\t" + fc;
		}
	}

}
