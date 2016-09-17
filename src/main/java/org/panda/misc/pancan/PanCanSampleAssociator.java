package org.panda.misc.pancan;

import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.TermCounter;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * @author Ozgun Babur
 */
public class PanCanSampleAssociator
{
	Map<String, String> sampleToType;
	Map<String, Integer> typeCounts;
	String[] matrixHeader;
	String matrixFile;

	public PanCanSampleAssociator() throws IOException
	{
		sampleToType = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/PanCan/SampleToCancerType.txt")).skip(1)
			.map(l -> l.split("\t"))
			.filter(t -> t[1].equals("Tumor"))
			.forEach(t -> sampleToType.put(t[0].substring(0, 15), t[2]));

		typeCounts = new HashMap<>();
		for (String type : sampleToType.values())
		{
			if (!typeCounts.containsKey(type)) typeCounts.put(type, 0);
			typeCounts.put(type, typeCounts.get(type) + 1);
		}

		matrixFile = "/home/babur/Documents/mutex/TCGA/PanCan/1/1/DataMatrix.txt";
		matrixHeader = Files.lines(Paths.get(matrixFile)).findFirst().get().split("\t");
	}

	public void printTypeDistribution(Collection<String> samples)
	{
		TermCounter tc = new TermCounter();
		tc.setGlobalCounts(typeCounts);
		for (String sample : samples)
		{
			tc.addTerm(getCancerTypeOf(sample));
		}
		tc.print();
	}

	public void writeSampleToType() throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(base + "SampleToType.txt"));
		writer.write("Sample\tType");
		Files.lines(Paths.get(base + "SampleToCancerType.txt")).skip(1).map(l -> l.split("\t"))
			.filter(t -> t[1].equals("Tumor")).forEach(t -> FileUtil.lnwrite(t[0].substring(0, 15) + "\t" + t[2], writer));
		writer.close();
	}

	public Map<String, Integer> getTypeCounts(Collection<String> samples)
	{
		Map<String, Integer> cnt = new HashMap<>();
		for (String sample : samples)
		{
			String type = sampleToType.get(sample);
			if (!cnt.containsKey(type)) cnt.put(type, 1);
			else cnt.put(type, cnt.get(type) + 1);
		}
		return cnt;
	}

	public String getCancerTypeOf(String sample)
	{
		return sampleToType.get(sample);
	}

	public Map<String, boolean[]> getTypePositions(String[] header)
	{
		assert !sampleToType.containsKey(header[0]) : "Samples should start from the second column";

		Map<String, Integer> sampleIndices = new HashMap<>();
		for (int i = 1; i < header.length; i++)
		{
			sampleIndices.put(header[i], i - 1);
		}

		Map<String, boolean[]> pos = new HashMap<>();

		for (int i = 1; i < header.length; i++)
		{
			String type = sampleToType.get(header[i]);
			if (!pos.containsKey(type)) pos.put(type, new boolean[header.length - 1]);
			pos.get(type)[sampleIndices.get(header[i])] = true;
		}
		return pos;
	}


	public Map<String, Set<String>> getAlteredSamples(String filename, Collection<String> genes) throws IOException
	{
		Map<String, Set<String>> geneToAlteredSamples = new HashMap<>();

		InputStream is = new FileInputStream(filename);
		if (filename.endsWith(".gz")) is = new GZIPInputStream(is);
		Scanner sc = new Scanner(new InputStreamReader(is));
		String[] header = sc.nextLine().split("\t");

		while (sc.hasNextLine())
		{
			String[] t = sc.nextLine().split("\t");
			if (genes != null && !genes.contains(t[0])) continue;
			Set<String> samples = new HashSet<>();
			for (int i = 1; i < t.length; i++)
			{
				if (!t[i].equals("0")) samples.add(header[i]);
			}
			geneToAlteredSamples.put(t[0], samples);
		}

//		Files.lines(Paths.get(filename != null ? filename : matrixFile)).skip(1).map(l -> l.split("\t"))
//			.filter(t -> genes == null || genes.contains(t[0])).forEach(t ->
//			{
//				Set<String> samples = new HashSet<>();
//				for (int i = 1; i < t.length; i++)
//				{
//					if (!t[i].equals("0")) samples.add(matrixHeader[i]);
//				}
//				geneToAlteredSamples.put(t[0], samples);
//			});

		return geneToAlteredSamples;
	}


	public static void main(String[] args) throws IOException
	{
//		List<String> genes = Arrays.asList("BTBD11", "RAI1");
//		PanCanSampleAssociator sa = new PanCanSampleAssociator();
//		Map<String, Set<String>> alteredSamples = sa.getAlteredSamples(sa.matrixFile, genes);
//		for (String gene : alteredSamples.keySet())
//		{
//			System.out.println("\ngene = " + gene);
//			Set<String> samples = alteredSamples.get(gene);
//			Map<String, Integer> counts = sa.getTypeCounts(samples);
//			counts.keySet().stream().sorted((t1, t2) -> counts.get(t2).compareTo(counts.get(t1))).forEach(t ->
//				System.out.println(t + "\t" + counts.get(t)));
//		}

		PanCanSampleAssociator pa = new PanCanSampleAssociator();
		pa.writeSampleToType();
	}
}
