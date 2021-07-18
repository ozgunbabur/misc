package org.panda.misc.proteomics;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class AddRegulatorsToPhosphoproteomicData
{
	public static final String PRIOR_FILE = "/Users/ozgun/Documents/data/causal-priors.txt";

	private static Map<String, Set<String>> posSU = new HashMap<>();
	private static Map<String, Set<String>> negSU = new HashMap<>();

	private static Map<String, Map<Integer, Set<String>>> posSS = new HashMap<>();
	private static Map<String, Map<Integer, Set<String>>> negSS = new HashMap<>();

	public static void main(String[] args) throws IOException
	{
		loadPriors();

		String inDir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/";
		String outDir = inDir + "regs-added/";
		Files.createDirectories(Paths.get(outDir));

		String[] files = new String[]{"human-phosphoproteome-ct-moderated.tsv",
			"human-phosphoproteome-evp-moderated.tsv",
			"mouse-phosphoproteome-pr443-evp-moderated.tsv",
			"mouse-phosphoproteome-pr607-ct-moderated.tsv",
			"mouse-phosphoproteome-pr607-evp-moderated.tsv"};

		for (String file : files)
		{
			addRegulators(inDir + file, outDir + file);
		}
	}

	private static void addRegulators(String sourceFile, String outFile) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		String[] header = Files.lines(Paths.get(sourceFile)).findFirst().get().split("\t");
		for (int i = 0; i < 4; i++)
		{
			writer.write(header[i] + "\t");
		}
		writer.write("Positive specific regulators\tNegative specific regulators\tPositive unspecific regulators\tNegative unspecific regulators");
		for (int i = 4; i < header.length; i++)
		{
			writer.write("\t" + header[i]);
		}

		Files.lines(Paths.get(sourceFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			List<String> genes = List.of(t[1].split(" "));
			List<String> siteStrs = List.of(t[2].split(" "));

			Set<String> posSpec = new HashSet<>();
			Set<String> negSpec = new HashSet<>();
			Set<String> posUnsp = new HashSet<>();
			Set<String> negUnsp = new HashSet<>();

			for (int i = 0; i < genes.size(); i++)
			{
				String gene = genes.get(i);

				if (posSU.containsKey(gene)) posUnsp.addAll(posSU.get(gene));
				if (negSU.containsKey(gene)) negUnsp.addAll(negSU.get(gene));

				for (String siteStr : siteStrs.get(i).split("\\|"))
				{
					Integer site = Integer.valueOf(siteStr.substring(1));

					if (posSS.containsKey(gene) && posSS.get(gene).containsKey(site))
					{
						posSpec.addAll(posSS.get(gene).get(site));
					}
					if (negSS.containsKey(gene) && negSS.get(gene).containsKey(site))
					{
						negSpec.addAll(negSS.get(gene).get(site));
					}
				}

				FileUtil.lnwrite(ArrayUtil.merge("\t", t[0], t[1], t[2], t[3]), writer);
				FileUtil.tab_write(CollectionUtil.merge(posSpec, " "), writer);
				FileUtil.tab_write(CollectionUtil.merge(negSpec, " "), writer);
				FileUtil.tab_write(CollectionUtil.merge(posUnsp, " "), writer);
				FileUtil.tab_write(CollectionUtil.merge(negUnsp, " "), writer);
				for (int j = 4; j < t.length; j++)
				{
					FileUtil.tab_write(t[j], writer);
				}
			}
		});

		writer.close();
	}

	private static void loadPriors() throws IOException
	{
		Files.lines(Paths.get(PRIOR_FILE)).map(l -> l.split("\t")).filter(t -> t[1].contains("phospho")).forEach(t ->
		{
			String source = t[0];
			String target = t[2];

			Map<String, Set<String>> su = t[1].startsWith("de") ? negSU : posSU;
			Map<String, Map<Integer, Set<String>>> ss = t[1].startsWith("de") ? negSS : posSS;

			if (!su.containsKey(target)) su.put(target, new HashSet<>());
			su.get(target).add(source);

			if (t.length > 4 && !t[4].isEmpty())
			{
				if (!ss.containsKey(target)) ss.put(target, new HashMap<>());

				for (String siteStr : t[4].split(";"))
				{
					Integer site = Integer.valueOf(siteStr.substring(1));

					if (!ss.get(target).containsKey(site)) ss.get(target).put(site, new HashSet<>());
					ss.get(target).get(site).add(source);
				}
			}
		});
	}
}
