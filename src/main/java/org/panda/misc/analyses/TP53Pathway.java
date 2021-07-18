package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.misc.MutexReader;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.StreamDirection;

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
//		readResults();
//		printBasalBRCATP53ComparisonCP();
//		cpResultSubgraph();
		printSignificantSize();
	}

	public static void generateMatrices() throws IOException
	{
//		Set<String> cases = Files.lines(Paths.get(BASE + "pancan_samples.txt")).skip(1).map(l -> l.split("\t")[0]).collect(Collectors.toSet());
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

	public static void printBasalBRCATP53ComparisonCP() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/TP53-CP/";
		Set<String> mutatedSamples = Files.lines(Paths.get(dir + "TP53-mutated-TCGA-BRCA-samples.tsv")).skip(1).map(l -> l.split("\t")[1].substring(0, 12)).collect(Collectors.toSet());
		Set<String> basalSamples = Files.lines(Paths.get(dir + "pancan_samples.txt")).skip(1).map(l -> l.split("\t")).filter(t -> t[2].equals("BRCA") && t[3].equals("Basal")).map(t -> t[0]).collect(Collectors.toSet());
		String[] header = Files.lines(Paths.get(dir + "CPTAC-TCGA-BRCA-data-77.txt")).findFirst().get().split("\t");

		for (String sample : header)
		{
			if (basalSamples.contains(sample))
			{
				if (mutatedSamples.contains(sample))
				{
					System.out.println("test-value-column = " + sample);
				}
				else
				{
					System.out.println("control-value-column = " + sample);
				}
			}
		}
	}

	public static void cpResultSubgraph() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/CPTAC-LUAD/mutations/TP53/";
		Set<String> goi = new HashSet<>(Arrays.asList(GOI.split("\n")));
		String subname = "goi-subset";
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-" + subname + ".sif", StreamDirection.DOWNSTREAM);
		Set<String> ids = CausalPathSubnetwork.getIDsAtTheNeighborhood(dir + "/results.txt", goi, StreamDirection.DOWNSTREAM);
		System.out.println("ids.size() = " + ids.size());
		CausalPathSubnetwork.writeSubsetFormat(dir + "/causative.format", dir + "/causative-" + subname + ".format", goi, ids);

	}

	public static void printSignificantSize() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/CPTAC-LUAD/mutations/TP53/";
		Set<String>[] sets = CausalPathSubnetwork.getSignificantGenes(dir + "significance-pvals.txt", 0.1);
		System.out.println("Act = " + sets[0].size());
		System.out.println("sets[0] = " + sets[0]);
		System.out.println("inh = " + sets[1].size() + "\n");
		System.out.println("sets[1] = " + sets[1]);


	}

	private static final String GOI =
		"CHEK1\n" +
		"ATR\n" +
		"ATM\n" +
		"CHEK2\n" +
		"TP53\n" +
		"MYC\n" +
		"MAX\n" +
		"CDK1\n" +
		"CDK2\n" +
		"RB1\n" +
		"RB2\n" +
		"RBL1\n" +
		"RBL2\n" +
		"MDM2\n" +
		"MDM4\n" +
		"PPM1D\n" +
		"CDKN2A\n" +
		"HUWE1";
}
