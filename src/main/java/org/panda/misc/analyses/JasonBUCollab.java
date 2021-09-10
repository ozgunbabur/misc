package org.panda.misc.analyses;

import org.panda.resource.MSigDB;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.GeneSetEnrichment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class JasonBUCollab
{
	public static final String DIR = "/Users/ozgun/Documents/Analyses/Jason/BU-collab/";

	public static void main(String[] args) throws IOException
	{
		Map<String, Double> signedP = readSignedP();
		Map<String, Double> pvals = removeSigns(signedP);

		double pThr = FDR.getPValueThreshold(pvals, null, 0.1);
		System.out.println("pThr = " + pThr);

		Set<String> highSkin = pvals.keySet().stream().filter(s -> pvals.get(s) <= pThr && signedP.get(s) < 0).collect(Collectors.toSet());
		Set<String> highVagina = pvals.keySet().stream().filter(s -> pvals.get(s) <= pThr && signedP.get(s) > 0).collect(Collectors.toSet());

		CollectionUtil.printVennSets(10, highSkin, highVagina);

		writeCPFile(signedP);

		doGSEA(highVagina, pvals.keySet(), "vagina-enrichment.tsv");
		doGSEA(highSkin, pvals.keySet(), "skin-enrichment.tsv");
		Set<String> allDif = new HashSet<>(highSkin);
		allDif.addAll(highVagina);
		doGSEA(allDif, pvals.keySet(), "all-enrichment.tsv");
	}

	static String getGN(String s)
	{
		int ind = s.indexOf(" GN=");
		if (ind > 0)
		{
			int end = s.indexOf(" ", ind + 4);
			return end > 0 ? s.substring(ind + 4, end) : s.substring(ind + 4);
		}
		return null;
	}

	static Map<String, Double> readTissue(String file)
	{
		Map<String, Double> map = new HashMap<>();

		String filename = DIR + file;
		String[] header = FileUtil.readHeader(filename);

		int descInd = ArrayUtil.indexOf(header, "Description");
		int pInd = ArrayUtil.indexOf(header, "p value(2-6)");

		FileUtil.linesTabbed(filename).skip(1).forEach(t ->
		{
			String sym = getGN(t[descInd]);
			if (sym == null) return;

			double p = Double.valueOf(t[pInd]);

			map.put(sym, Math.min(p, map.getOrDefault(sym, 1D)));
		});

		return map;
	}

	static Set<String> readAll()
	{
		return FileUtil.linesTabbed(DIR + "All.csv").skip(1).map(t -> getGN(t[2]))
			.filter(Objects::nonNull).collect(Collectors.toSet());
	}

	static void writeCPFile(Map<String, Double> signedP) throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(DIR + "vagina-vs-skin-causalpath/data.tsv");
		writer.write("ID\tSymbols\tSites\tEffects\tFeature\tSignedP");

		signedP.keySet().forEach(sym -> FileUtil.lnwrite(sym + "\t" + sym + "\t\t\tG\t" + signedP.get(sym), writer));

		writer.close();
	}

	static void doGSEA(Set<String> selected, Set<String> background, String outFile) throws IOException
	{
//		Map<String, Set<String>> hallmark = MSigDB.get().getSetsNameFiltered(name -> name.startsWith("HALLMARK_"));
		Map<String, Set<String>> reactome = MSigDB.get().getSetsNameFiltered(name -> name.startsWith("REACTOME_"));

		Map<String, Double> pvals = GeneSetEnrichment.calculateEnrichment(selected, background, 5, 200, reactome);
		GeneSetEnrichment.writeEnrichment(selected, reactome, pvals, DIR + outFile);
	}

	static Map<String, Double> readSignedP()
	{
		String filename = DIR + "SkinVsVagTake2.csv";
		String[] header = FileUtil.readHeader(filename);

		int descInd = ArrayUtil.indexOf(header, "Description");
		int pInd = ArrayUtil.indexOf(header, "p value(S/V)");
		int dInd = ArrayUtil.indexOf(header, "expressed more");

		return FileUtil.linesTabbed(filename).skip(1).peek(t -> t[descInd] = getGN(t[descInd]))
			.filter(t -> t[descInd] != null).collect(Collectors.toMap(t -> t[descInd], t -> Double.valueOf(t[pInd]) *
				(t[dInd].equals("v") ? 1 : -1), (v1, v2) -> (Math.abs(v1) < Math.abs(v2) ? v1 : v2)));
	}

	static Map<String, Double> removeSigns(Map<String, Double> map)
	{
		return map.keySet().stream().collect(Collectors.toMap(k -> k, k -> Math.abs(map.get(k))));
	}
}
