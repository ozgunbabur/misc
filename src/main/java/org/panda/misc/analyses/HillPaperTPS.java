package org.panda.misc.analyses;

import org.panda.resource.UniProtSequence;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class HillPaperTPS
{
	public static final String INPUT_BASE = "/Users/ozgun/Documents/Analyses/CausalPath-data/CL-Ligand-Drug/DataS1/converted/";
	public static final String OUTPUT_BASE = "/Users/ozgun/Documents/Analyses/CausalPath-data/CL-Ligand-Drug-TPS/";
	public static final String CP_BASE = "/Users/ozgun/Documents/Analyses/CausalPath-data/CL-Ligand-Drug/";


	public static final Map<String, String> sourcesMap = new HashMap<>();
	static
	{
		sourcesMap.put("EGF", "EGF_HUMAN");
		sourcesMap.put("FGF1", "FGF1_HUMAN");
		sourcesMap.put("HGF", "HGF_HUMAN");
		sourcesMap.put("IGF1", "IGF1_HUMAN");
		sourcesMap.put("Insulin", "INSR_HUMAN");
		sourcesMap.put("NRG1", "NRG1_HUMAN");
	}

	public static void main(String[] args) throws IOException
	{
//		prepareTimeSeriesInput();
//		prepareOtherInputs();
//		writeAllExpectedChanges(OUTPUT_BASE, OUTPUT_BASE + "expected.txt");
		testHypotheses(OUTPUT_BASE + "expected.txt", CP_BASE + "DataS1/converted/");
	}

	private static void prepareTimeSeriesInput() throws IOException
	{
		for (File cFile : new File(INPUT_BASE).listFiles())
		{
			if (cFile.getName().startsWith(".")) continue;

			String cell = cFile.getName();
			cell = cell.substring(0, cell.indexOf("."));

			String[] header = Files.lines(Paths.get(cFile.getPath())).findFirst().get().split("\t");

			Map<String, Map<Integer, String>> map = new HashMap<>();

			for (int i = 4; i < header.length; i++)
			{
				if (header[i].startsWith("-") || header[i].startsWith("PBS") || !header[i].contains("-DMSO-")) continue;

				String[] t = header[i].split("-");

				String ligand = t[0];
				String time = t[2];
				int min = timeInMins(time);

				if (!map.containsKey(ligand)) map.put(ligand, new HashMap<>());
				map.get(ligand).put(min, header[i]);
			}

			Map<Integer, String> ctrlMap = map.get("Serum");

			for (String ligand : map.keySet())
			{
				if (ligand.equals("Serum")) continue;

				Map<Integer, String> ligandMap = map.get(ligand);

				List<List<Integer>> groups = partitionInTriplets(ligandMap.keySet());

				List<List<String>> samples = new ArrayList<>();
				List<String> initial = new ArrayList<>();
				for (Integer min : groups.get(0))
				{
					initial.add(ctrlMap.get(min));
				}
				samples.add(initial);

				for (List<Integer> group : groups)
				{
					List<String> sampleGroup = new ArrayList<>();
					for (Integer min : group)
					{
						sampleGroup.add(ligandMap.get(min));
					}
					samples.add(sampleGroup);
				}

				String outBase = OUTPUT_BASE + cell + "/" + ligand + "/timeseries/";
				Files.createDirectories(Paths.get(outBase));

				BufferedWriter medianWriter = Files.newBufferedWriter(Paths.get(outBase + "median-time-series.tsv"));
				BufferedWriter firstWriter = Files.newBufferedWriter(Paths.get(outBase + "p-values-first.tsv"));
				BufferedWriter prevWriter = Files.newBufferedWriter(Paths.get(outBase + "p-values-prev.tsv"));
				BufferedWriter mapWriter = Files.newBufferedWriter(Paths.get(outBase + "peptide-mapping.tsv"));

				mapWriter.write("peptide\tprotein(s)");
				String medianHeader = "peptide";
				for (int i = 0; i < samples.size(); i++)
				{
					medianHeader += "\t" + CollectionUtil.merge(samples.get(i), "|");
				}
				medianWriter.write(medianHeader);
				firstWriter.write("#peptide");
				prevWriter.write("#peptide");
				String[] tt = medianHeader.split("\t");
				for (int i = 2; i < tt.length; i++)
				{
					firstWriter.write("\t" + tt[i] + " vs " + tt[1]);
					prevWriter.write("\t" + tt[i] + " vs " + tt[i-1]);
				}

				Files.lines(Paths.get(cFile.getPath())).skip(1).map(l -> l.split("\t")).forEach(t ->
				{
					String id = t[0];
					String symbols = t[1];
					String names = getUPNamesFromHGNC(symbols);
					FileUtil.lnwrite(id + "\t" + names, mapWriter);

					List<double[]> valsList = new ArrayList<>();

					for (int i = 0; i < samples.size(); i++)
					{
						List<String> group = samples.get(i);
						double[] arr = new double[group.size()];
						valsList.add(arr);
						for (int j = 0; j < group.size(); j++)
						{
							String sample = group.get(j);
							double val = Double.valueOf(t[ArrayUtil.indexOf(header, sample)]);
							arr[j] = val;
						}
					}

					FileUtil.lnwrite(id, medianWriter);
					FileUtil.lnwrite(id, firstWriter);
					FileUtil.lnwrite(id, prevWriter);

					for (int i = 0; i < valsList.size(); i++)
					{
						double median = Summary.median(valsList.get(i));
						FileUtil.tab_write(median, medianWriter);

						if (i > 0)
						{
							Tuple testFirst = TTest.test(valsList.get(0), valsList.get(i));
							Tuple testPrev = TTest.test(valsList.get(i-1), valsList.get(i));

							FileUtil.tab_write(testFirst.p, firstWriter);
							FileUtil.tab_write(testPrev.p, prevWriter);
						}
					}
				});

				medianWriter.close();
				firstWriter.close();
				prevWriter.close();
				mapWriter.close();
			}

		}
	}

	private static void prepareOtherInputs() throws IOException
	{
		for (File cellDir : new File(OUTPUT_BASE).listFiles())
		{
			if (cellDir.getName().startsWith(".") || !cellDir.isDirectory()) continue;

			for (File ligDir : cellDir.listFiles())
			{
				if (ligDir.getName().startsWith(".") || !ligDir.isDirectory()) continue;

				String ligand = ligDir.getName();

				String pcsfDir = ligDir.getPath() + "/pcsf/";
				Files.createDirectories(Paths.get(pcsfDir));

				BufferedWriter writer = Files.newBufferedWriter(Paths.get(pcsfDir + "sources.txt"));
				writer.write(sourcesMap.get(ligand));
				writer.close();

				Files.createDirectories(Paths.get(ligDir.getPath() + "/networks"));

				writer = Files.newBufferedWriter(Paths.get(ligDir.getPath() + "/generate_prizes.sh"));
				writer.write(getGeneratePrizesScript());
				writer.close();

				writer = Files.newBufferedWriter(Paths.get(ligDir.getPath() + "/generate_PCSF_network.sh"));
				writer.write(getPCSFScript());
				writer.close();

				writer = Files.newBufferedWriter(Paths.get(ligDir.getPath() + "/runTPS.sh"));
				writer.write(getRunTPSScript(ligand, ligDir.getPath()));
				writer.close();

				writer = Files.newBufferedWriter(Paths.get(ligDir.getPath() + "/run.sh"));
				writer.write(getRunScript());
				writer.close();

				Runtime.getRuntime().exec("chmod 755 " + ligDir.getPath() + "/*.sh");
			}
		}
	}

	private static String getGeneratePrizesScript()
	{
		return "#!/bin/bash\n" +
			"# A wrapper script to generate protein prizes from the TPS peptide-protein\n" +
			"# map, first scores file, and previous scores file.\n" +
			"\n" +
			"python2 /home/ozgun/Documents/Code/tps/pcsf/generate_prizes.py \\\n" +
			"    --firstfile=timeseries/p-values-first.tsv \\\n" +
			"\t--prevfile=timeseries/p-values-prev.tsv \\\n" +
			"\t--mapfile=timeseries/peptide-mapping.tsv \\\n" +
			"\t--outfile=pcsf/prizes.txt\n";
	}

	private static String getPCSFScript()
	{
		return "#!/bin/bash\n" +
			"# Submit PCSF runs with different seeds for random edge noise\n" +
			"# Sets up a new configuration file if needed and uses\n" +
			"# environment variables to pass other arguments to forest.py\n" +
			"\n" +
			"# Set the output directory for PCSF and HTCondor output\n" +
			"export outpath=pcsf-results\n" +
			"mkdir -p $outpath\n" +
			"\n" +
			"# Set the code paths for the Omics Integrator and msgsteiner dependencies\n" +
			"## This is the directory that contains scripts/forest.py\n" +
			"export oipath=/home/ozgun/Documents/OmicsIntegrator-0.3.1\n" +
			"## This is the path to the msgsteiner executable, including the executable name\n" +
			"export msgsteinerpath=/home/ozgun/Documents/msgsteiner-1.3/msgsteiner\n" +
			"\n" +
			"# Parameters set during the parameter sweep\n" +
			"b=0.55\n" +
			"m=0.008\n" +
			"w=0.1\n" +
			"\n" +
			"# Prize filename prefix (assume a .txt extension follows, e.g. firstprev.txt)\n" +
			"# It will be used to create an output prefix for the Steiner forest networks\n" +
			"## This example uses the prizes derived from the statistical significance of each\n" +
			"## protein's phosphorlyation intensity after EGF stimulation compared to its\n" +
			"## phosphorylation at the previous and first time points\n" +
			"export prizetype=prizes\n" +
			"# The path to the prize file above\n" +
			"export prizepath=pcsf\n" +
			"# The PPI network, including the path\n" +
			"## This example uses a combination of PhosphoSitePlus and iRefIndex interactions\n" +
			"## with UniProt entry name identifiers\n" +
			"export edgefile=/home/ozgun/Documents/Code/tps/data/networks/phosphosite-irefindex13.0-uniprot.txt\n" +
			"# A file listing the protein names that should be treated as source nodes,\n" +
			"# including the path\n" +
			"export sources=pcsf/sources.txt\n" +
			"\n" +
			"# The following three parameters can typically be left at these default values\n" +
			"# Depth from root of tree\n" +
			"D=10\n" +
			"# Convergence parameter\n" +
			"g=1e-3\n" +
			"# Noise to compute a family of solutions\n" +
			"r=0.01\n" +
			"\n" +
			"# Create the configuration file, removing an older copy of the file if it exists\n" +
			"mkdir -p ${outpath}/conf\n" +
			"filename=${outpath}/conf/conf_w${w}_b${b}_D${D}_m${m}_r${r}_g${g}.txt\n" +
			"rm -f $filename\n" +
			"touch $filename\n" +
			"printf \"w = ${w}\\n\" >> $filename\n" +
			"printf \"b = ${b}\\n\" >> $filename\n" +
			"printf \"D = ${D}\\n\" >> $filename\n" +
			"printf \"mu = ${m}\\n\" >> $filename\n" +
			"printf \"r = ${r}\\n\" >> $filename\n" +
			"printf \"g = ${g}\\n\" >> $filename\n" +
			"\n" +
			"# Set the remaining environment variables\n" +
			"export conf=$filename\n" +
			"export beta=$b\n" +
			"export mu=$m\n" +
			"export omega=$w\n" +
			"\n" +
			"# Use different seeds for each run, which will control the random edge noise\n" +
			"# Create a family of 100 forests\n" +
			"for s in $(seq 1 100)\n" +
			"do\n" +
			"\t# Set the seed\n" +
			"\texport seed=$s\n" +
			"\n" +
			"\t# Submit the job to HTCondor with the configuration file, params, and seed\n" +
			"\t# Could replace this with a submission to a different queueing system\n" +
			"\t# (e.g. qsub instead of condor_submit) or directly call run_PCSF.sh to\n" +
			"\t# run locally.\n" +
			"\t/home/ozgun/Documents/Code/tps/pcsf/run_PCSF.sh\n" +
			"done\n" +
			"\n" +
			"# Generate a wrapper script to summarize the family of forests\n" +
			"# This must be run after all PCSF runs terminate\n" +
			"# HTCondor can manage these dependencies with the Directed Acyclic Graph Manager\n" +
			"# but this strategy generalizes to other setups\n" +
			"# Set the name of the summarization script, which overwrites an existing\n" +
			"# file with the same name\n" +
			"# The script assumes that the summarization Python code resides in the same\n" +
			"# directory\n" +
			"sumscript=summarize_forests.sh\n" +
			"# A filename pattern used to collect all of the forest output files\n" +
			"pattern=${prizetype}_beta${beta}_mu${mu}_omega${omega}\n" +
			"rm -f $sumscript\n" +
			"touch $sumscript\n" +
			"printf \"#!/bin/bash\\n\" >> $sumscript\n" +
			"printf \"#Summarize a family of Steiner forests\\n\" >> $sumscript\n" +
			"printf \"python2 /home/ozgun/Documents/Code/tps/pcsf/summarize_sif.py --indir ${outpath} --pattern ${pattern}*optimalForest.sif --prizefile ${prizepath}/${prizetype}.txt --outfile ${outpath}/${pattern}_summary\\n\" >> $sumscript\n" +
			"chmod u+x $sumscript\n" +
			"\n";
	}

	private static String getRunTPSScript(String ligand, String path)
	{
		path = path.replace("/Users/", "/home/");
		return "cd /home/ozgun/Documents/Code/tps\n" +
			"rm output.sif\n" +
			"rm output.tsv\n" +
			"rm temporal-interpretation.tsv\n" +
			"rm activity-windows.tsv\n\n" +
			"./scripts/run \\\n" +
			"   --network " + path + "/networks/input-network.tsv \\\n" +
			"   --timeseries " + path + "/timeseries/median-time-series.tsv \\\n" +
			"   --firstscores " + path + "/timeseries/p-values-first.tsv \\\n" +
			"   --prevscores " + path + "/timeseries/p-values-prev.tsv \\\n" +
			"   --peptidemap " + path + "/timeseries/peptide-mapping.tsv \\\n" +
			"   --source " + sourcesMap.get(ligand) + " \\\n" +
			"   --threshold 0.01\n\n" +
			"cp output.sif " + path + "/.\n" +
			"cp output.tsv " + path + "/.\n" +
			"cp temporal-interpretation.tsv " + path + "/.\n" +
			"cp activity-windows.tsv " + path + "/.\n";
	}

	private static String getRunScript()
	{
		return "#!/bin/bash\n" +
			"sh generate_prizes.sh\n" +
			"sh generate_PCSF_network.sh\n" +
			"sh summarize_forests.sh\n" +
			"cp pcsf-results/prizes_beta0.55_mu0.008_omega0.1_summary_union.tsv networks/input-network.tsv\n" +
			"sh runTPS.sh\n";
	}

	private static int timeInMins(String time)
	{
		if (time.endsWith("min")) return Integer.valueOf(time.substring(0, time.length() - 3));
		if (time.endsWith("hr")) return Integer.valueOf(time.substring(0, time.length() - 2)) * 60;
		throw new RuntimeException("Illegal time string: " + time);
	}

	private static List<List<Integer>> partitionInTriplets(Set<Integer> timepoints)
	{
		int numGroups = timepoints.size() / 3;

		int irreg = timepoints.size() % 3;

		List<Integer> groupSizes = new ArrayList<>();
		for (int i = 0; i < numGroups; i++)
		{
			int size = 3;
			if (irreg > i) size++;
			groupSizes.add(size);
		}

		List<Integer> sorted = new ArrayList<>(timepoints);
		sorted.sort(Integer::compareTo);

		List<List<Integer>> groups = new ArrayList<>();

		int nextInd = 0;
		for (Integer size : groupSizes)
		{
			groups.add(sorted.subList(nextInd, nextInd + size));
			nextInd += size;
		}

		return groups;
	}

	private static String getUPNamesFromHGNC(String syms)
	{
		List<String> upNames = new ArrayList<>();
		for (String sym : syms.split(" "))
		{
			String name = UniProtSequence.get().getNameOfSymbol(sym, "9606");
			if (name != null) upNames.add(name);
		}
		return CollectionUtil.merge(upNames, "|");
	}

	private static Map<String, Set<String>> readDrugToTargets() throws IOException
	{
		return Files.lines(Paths.get(CP_BASE + "drug-targets.txt")).skip(1)
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0],
				t -> new HashSet<>(Arrays.asList(t[1].split(" ")))));
	}

	private static void writeAllExpectedChanges(String baseDir, String outFile) throws IOException
	{
		Map<String, Map<String, Map<String, Map<String, Map<Integer, String>>>>> map = getAllExpectedChanges(baseDir);
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		writer.write("Cell-line\tLigand\tDrug\tPepID\tExpected-direction\tTime-slots");

		for (String cell : map.keySet())
		{
			for (String ligand : map.get(cell).keySet())
			{
				for (String drug : map.get(cell).get(ligand).keySet())
				{
					for (String pepID : map.get(cell).get(ligand).get(drug).keySet())
					{
						map.get(cell).get(ligand).get(drug).get(pepID).forEach((ch, times) ->
							FileUtil.lnwrite(ArrayUtil.merge("\t", cell, ligand, drug, pepID,
								ch.toString(), times), writer));
					}
				}
			}
		}

		writer.close();
	}

	private static Map<String, Map<String, Map<String, Map<String, Map<Integer, String>>>>> getAllExpectedChanges(String baseDir) throws IOException
	{
		// replace gene symbols with UP names
		Map<String, Set<String>> drugToTargets = readDrugToTargets();
		for (String drug : new HashSet<>(drugToTargets.keySet()))
		{
			Set<String> symSet = drugToTargets.get(drug);
			Set<String> nameSet = symSet.stream().map(t -> UniProtSequence.get().getNameOfSymbol(t, "9606")).collect(Collectors.toSet());
			drugToTargets.put(drug, nameSet);
		}

		// cell-line --> ligand --> drug --> dwnID --> expected-change --> time-slot
		Map<String, Map<String, Map<String, Map<String, Map<Integer, String>>>>> map = new HashMap<>();

		for (File cellDir : new File(baseDir).listFiles())
		{
			if (cellDir.getName().startsWith(".") || !cellDir.isDirectory()) continue;

			map.put(cellDir.getName(), new HashMap<>());

			for (File ligandDir : cellDir.listFiles())
			{
				String resultFile = ligandDir.getPath() + "/output.sif";

				if (Files.exists(Paths.get(resultFile)))
				{
					map.get(cellDir.getName()).put(ligandDir.getName(), determineExpectationsAfterInhibitions(
						ligandDir.getPath(), drugToTargets));
				}
			}
		}

		return map;
	}

	private static Map<String, Map<String, Map<Integer, String>>> determineExpectationsAfterInhibitions(String resultDir, Map<String, Set<String>> drugToTargets) throws IOException
	{
		// drug --> downstreamID --> expected-change --> time-slots
		Map<String, Map<String, Map<Integer, String>>> map = new HashMap<>();

		Map<String, Set<String>> tgtToDrugs = reverse(drugToTargets);

		String[] header = Files.lines(Paths.get(resultDir + "/activity-windows.tsv")).findFirst().get().split("\t");

		// find the activation time slots of drug targets

		Map<String, Integer> firstActivatedIndexMap = new HashMap<>();

		Files.lines(Paths.get(resultDir + "/activity-windows.tsv")).skip(1).map(l -> l.split("\t"))
			.filter(t -> tgtToDrugs.containsKey(t[0].split("#")[0])).forEach(t ->
		{
			for (int i = 1; i < t.length; i++)
			{
				if (t[i].equals("activation"))
				{
					firstActivatedIndexMap.put(t[0].split("#")[0], i);
					break;
				}
			}
		});

		// find the most significant change of the downstream genes

		Map<String, Set<String>> tgtToDwn = new HashMap<>();

		Files.lines(Paths.get(resultDir + "/output.sif")).map(l -> l.split("\t")).filter(t -> !t[2].equals("U")).forEach(t ->
		{
			if (firstActivatedIndexMap.keySet().contains(t[0]))
			{
				if (!tgtToDwn.containsKey(t[0])) tgtToDwn.put(t[0], new HashSet<>());
				tgtToDwn.get(t[0]).add(t[2]);
			}
		});

		Map<String, Set<String>> dwnToTrg = reverse(tgtToDwn);

		Map<String, Set<String>> dwnToID = new HashMap<>();
		Files.lines(Paths.get(resultDir + "/timeseries/peptide-mapping.tsv")).skip(1)
			.map(l -> l.split("\t")).filter(t -> dwnToTrg.containsKey(t[1])).forEach(t ->
		{
			if (!dwnToID.containsKey(t[1])) dwnToID.put(t[1], new HashSet<>());
			dwnToID.get(t[1]).add(t[0]);
		});

		for (String dwn : dwnToID.keySet())
		{
			Set<String> trgs = dwnToTrg.get(dwn);
			for (String trg : trgs)
			{
				String[] sigEff = getMostSignificantEffect(resultDir + "/timeseries/p-values-first.tsv", dwnToID.get(dwn), firstActivatedIndexMap.get(trg));

				String id = sigEff[0];
				Integer time = Integer.valueOf(sigEff[1]);

				for (String drug : tgtToDrugs.get(trg))
				{
					Integer direction = getTheReverseChangeDirection(resultDir + "/timeseries/median-time-series.tsv", id, time);

					// put it in map
					if (!map.containsKey(drug)) map.put(drug, new HashMap<>());
					if (!map.get(drug).containsKey(id)) map.get(drug).put(id, new HashMap<>());
					map.get(drug).get(id).put(direction, header[time]);
				}
			}
		}

		return map;
	}

	private static Map<String, Set<String>> reverse(Map<String, Set<String>> map)
	{
		Map<String, Set<String>> rev = new HashMap<>();
		map.forEach((key, set) -> set.forEach(s ->
		{
			if (!rev.containsKey(s)) rev.put(s, new HashSet<>());
			rev.get(s).add(key);
		}));
		return rev;
	}

	/**
	 * Returns {ID, timeIndex} tuple that indicated the most significant change for the downstream protein. We will test
	 * the reversal of this change using inhibition experiments.
	 */
	private static String[] getMostSignificantEffect(String firstFile, Set<String> ids, int firstActivatedIndex) throws IOException
	{
		Map<String[], Double> pMap = new HashMap<>();

		Files.lines(Paths.get(firstFile)).skip(1).map(l -> l.split("\t"))
			.filter(t -> ids.contains(t[0])).forEach(t ->
		{
			for (int i = firstActivatedIndex; i < t.length; i++)
			{
				pMap.put(new String[]{t[0], "" + i}, Double.valueOf(t[i]));
			}
		});

		double best = pMap.values().stream().min(Double::compareTo).get();

		return pMap.keySet().stream().filter(k -> pMap.get(k).equals(best)).findFirst().get();
	}

	private static Integer getTheReverseChangeDirection(String medianFile, String id, int time) throws IOException
	{
		String[] token = Files.lines(Paths.get(medianFile)).skip(1).map(l -> l.split("\t"))
			.filter(t -> t[0].equals(id)).findFirst().get();

		return Double.valueOf(token[time + 1]) > Double.valueOf(token[1]) ? -1 : 1;
	}


	private static void testHypotheses(String expectationFile, String dataDir) throws IOException
	{
		Map<Tuple, Double> pvals = new HashMap<>();

		Files.lines(Paths.get(expectationFile)).skip(1).filter(l -> l.contains("_")).map(l -> l.split("\t")).forEach(hy ->
		{
			String cell = hy[0];
			String ligand = hy[1];
			String inh = hy[2];
			String id = hy[3];
			Integer exp = Integer.valueOf(hy[4]);
			String sampCue = hy[5];

			String dataFile = dataDir + cell + ".txt";
			String[] header = FileUtil.lines(dataFile).findFirst().get().split("\t");

			double[] values = ArrayUtil.readToDoubleArrayPreserveMissing(FileUtil.lines(dataFile)
				.filter(l -> l.startsWith(id + "\t")).findFirst().get().split("\t"));

			String ctrlPrefix = ligand + "-DMSO-";
			String testPrefix = ligand + "-" + inh + "-";
			Set<String> validSuffixes = getValidSuffixes(sampCue);

			Set<String> ctrlNames = getColNames(header, ctrlPrefix, validSuffixes);
			Set<String> testNames = getColNames(header, testPrefix, validSuffixes);

			Map<String, String> pairing = getPairing(ctrlNames, testNames);

			if (pairing.size() < 3) return;

			double[][] ct = getSubsets(values, header, pairing);

			Tuple tuple = TTest.testPaired(ct[0], ct[1]);

			if (tuple.isNaN()) return;

			tuple.v = getMeanChangePaired(ct[0], ct[1]);
			tuple.v /= Summary.stdev(values);
			tuple.v *= exp;

			if (-Math.log(tuple.p) > 8 && tuple.v > 0 && tuple.v < 0.2)
			{
				System.out.println();
			}

			pvals.put(tuple, tuple.p);

			System.out.println(ArrayUtil.getString("\t", cell, ligand, inh, id, exp, sampCue) + "\t" + tuple.v + "\t" + (-Math.log(tuple.p)));
		});

		List<Tuple> select = FDR.select(pvals, null, 0.1);
		System.out.println("select.size() = " + select.size());
		System.out.println("Agree: " + select.stream().filter(t -> t.v > 0).count());
		System.out.println("Disagree: " + select.stream().filter(t -> t.v < 0).count());

		double pValueThreshold = FDR.getPValueThreshold(pvals, null, 0.1);
		double logP = Math.log(pValueThreshold);
		System.out.println("pValueThreshold = " + pValueThreshold);
		System.out.println("-logP = " + -logP);

		Histogram his = new Histogram(1);
		his.setBorderAtZero(true);
		pvals.keySet().forEach(t -> his.count(t.v));
		his.print();
	}

	private static Set<String> getValidSuffixes(String s)
	{
		return Arrays.stream(s.split("\\|"))
			.map(sample -> sample.substring(sample.lastIndexOf("-"))).collect(Collectors.toSet());
	}

	private static Set<String> getColNames(String[] header, String prefix, Set<String> validSuffixes)
	{
		Set<String> names = new HashSet<>();
		for (int i = 0; i < header.length; i++)
		{
			int ii = i;
			if (header[i].startsWith(prefix) && validSuffixes.stream().anyMatch(s -> header[ii].endsWith(s))) names.add(header[i]);
		}
		return names;
	}

	private static Map<String, String> getPairing(Set<String> controls, Set<String> tests)
	{
		Map<String, String> map = new HashMap<>();

		for (String control : controls)
		{
			String suffix = control.substring(control.lastIndexOf("-") + 1);
			String chosen = null;

			for (String test : tests)
			{
				if (test.substring(test.lastIndexOf("-") + 1).equals(suffix))
				{
					chosen = test;
					break;
				}
			}

			if (chosen != null) map.put(control, chosen);
		}
		return map;
	}

	private static double[][] getSubsets(double[] vals, String[] header, Map<String, String> pairing)
	{
		double[] ctrl = new double[pairing.size()];
		double[] test = new double[pairing.size()];

		int j = 0;

		for (String contName : pairing.keySet())
		{
			ctrl[j]   = vals[ArrayUtil.indexOf(header, contName)];
			test[j++] = vals[ArrayUtil.indexOf(header, pairing.get(contName))];
		}

		return new double[][]{ctrl, test};
	}

	private static double getMeanChangePaired(double[] ctrl, double[] test)
	{
		double[] m = new double[test.length];
		for (int i = 0; i < m.length; i++)
		{
			m[i] = test[i] - ctrl[i];
		}
		return Summary.mean(m);
	}
}
