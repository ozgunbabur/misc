package org.panda.misc.analyses;

import org.panda.misc.MutexReader;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class TP53Pathway
{
	public static final String DIR = "/Users/ozgun/Documents/Analyses/TP53-mutex/";

	public static void main(String[] args) throws IOException
	{
//		generateMatrices();
		readResults();
	}

	public static void generateMatrices() throws IOException
	{
//		Set<String> cases = Files.lines(Paths.get(DIR + "pancan_samples.txt")).skip(1).map(l -> l.split("\t")[0]).collect(Collectors.toSet());
//		System.out.println("cases.size() = " + cases.size());

		Map<String, Set<String>> studyMap = new HashMap<>();
		Files.lines(Paths.get(DIR + "pancan_samples.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (!studyMap.containsKey(t[2])) studyMap.put(t[2], new HashSet<>());
			studyMap.get(t[2]).add(t[0]);
		});

		System.out.println("studyMap.size() = " + studyMap.size());


		Map<String, Set<String>> caseToGenes = new HashMap<>();

		Files.lines(Paths.get(DIR + "TP53_AssocGenehotspot_input.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (!caseToGenes.containsKey(t[1])) caseToGenes.put(t[1], new HashSet<>());
			caseToGenes.get(t[1]).add(t[0]);
		});

		Files.lines(Paths.get(DIR + "TP53count_mut_tcgaHs.txt")).skip(1).map(l -> l.split("\t")[0]).forEach(c ->
		{
			if (!caseToGenes.containsKey(c)) caseToGenes.put(c, new HashSet<>());
			caseToGenes.get(c).add("TP53");
		});


		for (String study : studyMap.keySet())
		{
			String dir = DIR + study;
			Files.createDirectories(Paths.get(dir));

			BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "/matrix.txt"));
			writer.write("Gene");

			Set<String> caseNames = studyMap.get(study);
			Set<String> genes = caseToGenes.keySet().stream().filter(caseNames::contains).map(caseToGenes::get).flatMap(Collection::stream).collect(Collectors.toSet());

			caseNames.stream().sorted().forEach(c -> FileUtil.tab_write(c, writer));

			genes.stream().sorted().forEach(g ->
			{
				FileUtil.lnwrite(g, writer);

				caseNames.stream().sorted().forEach(c -> FileUtil.tab_write(caseToGenes.getOrDefault(c, Collections.emptySet()).contains(g) ? "1" : "0", writer));
			});

			writer.close();
		}
	}

	public static void readResults()
	{
		HashSet<MutexReader.Group> results = new HashSet<>();
		MutexReader.readMutexResultsRecursive(DIR, results);

		results.stream().sorted(Comparator.comparingDouble(o -> o.score)).filter(g -> g.score < 1).forEach(g -> System.out.println(g.toString().replace(DIR, "")));
	}
}
