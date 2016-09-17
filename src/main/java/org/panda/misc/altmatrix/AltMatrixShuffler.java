package org.panda.misc.altmatrix;

import org.panda.utility.ArrayUtil;
import org.panda.utility.Progress;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * For randomizing an alteration matrix while preserving sample weights.
 * @author Ozgun Babur
 */
public class AltMatrixShuffler
{
	private int[][] matrix;
	private String header;
	private List<String> genes;

	private Map<String, int[]> subgroupIndex;
	private Map<Integer, String> indexToSubgroup;

	private Random rand = new Random();

	List<Integer> geneIndices;
	int currentRandomGeneIndex = -1;

	private void loadMatrix(String filename) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));
		header = sc.nextLine();
		genes = new ArrayList<>();
		List<int[]> vals = new ArrayList<>();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			genes.add(token[0]);
			int[] v = new int[token.length - 1];
			for (int i = 1; i < token.length; i++)
			{
				v[i - 1] = Integer.parseInt(token[i]);
			}
			vals.add(v);
		}
		matrix = vals.toArray(new int[vals.size()][]);
	}

	private void loadSubtypes(String filename) throws IOException
	{
		Map<String, String> sampleToType = Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> t[1], (s, s2) -> s));

		String[] colNames = header.substring(header.indexOf("\t") + 1).split("\t");

		Map<String, List<Integer>> tissIndices = sampleToType.values().stream().distinct()
			.collect(Collectors.toMap(t -> t, t -> new ArrayList<>()));

		for (int i = 0; i < colNames.length; i++)
		{
			tissIndices.get(sampleToType.get(colNames[i])).add(i);
		}

		subgroupIndex = new HashMap<>();
		tissIndices.keySet().forEach(t -> subgroupIndex.put(t, ArrayUtil.convertToBasicIntArray(tissIndices.get(t))));

		indexToSubgroup = new HashMap<>();
		subgroupIndex.keySet().forEach(t -> Arrays.stream(subgroupIndex.get(t)).forEach(i -> indexToSubgroup.put(i, t)));
	}

	private void shuffle()
	{
		Progress p = new Progress(matrix.length, "Shuffling");
		for (int g1 = 0; g1 < matrix.length; g1++)
		{
			for (int s1 = 0; s1 < matrix[g1].length; s1++)
			{
				int s2 = randSample(s1);

				if (contrasts(matrix[g1][s1], matrix[g1][s2]))
				{
					int g2 = findSuitableGene(s1, s2, matrix[g1][s1], matrix[g1][s2], g1);
					if (g2 >= 0)
					{
						swap(g1, s1, s2);
						swap(g2, s1, s2);
					}
				}
			}
			p.tick();
		}
	}

	private void swap(int g, int s1, int s2)
	{
		int temp = matrix[g][s1];
		matrix[g][s1] = matrix[g][s2];
		matrix[g][s2] = temp;
	}

	private int findSuitableGene(int s1, int s2, int alt1, int alt2, int currentGene)
	{
		for (int i = 0; i < matrix.length; i++)
		{
			int g = nextRandomGene();

			if (g != currentGene &&
				contrasts(matrix[g][s1], alt1) && contrasts(matrix[g][s2], alt2))
			{
				return g;
			}
		}
		return -1;
	}

	private int nextRandomGene()
	{
		if (geneIndices == null)
		{
			geneIndices = IntStream.range(0, matrix.length).boxed().collect(Collectors.toList());
			Collections.shuffle(geneIndices);
		}
		currentRandomGeneIndex++;
		if (currentRandomGeneIndex == geneIndices.size()) currentRandomGeneIndex = 0;

		return geneIndices.get(currentRandomGeneIndex);
	}

	private boolean contrasts(int i, int j)
	{
		return i + j > 0 && i * j == 0;
	}

	private int randSample(int current)
	{
		if (subgroupIndex == null) return rand(matrix[0].length, current);
		else
		{
			String sub = indexToSubgroup.get(current);
			int[] indices = subgroupIndex.get(sub);

			if (indices.length == 1) return current;

			int x;
			do
			{
				x = rand.nextInt(indices.length);
			}
			while (indices[x] == current);

			return indices[x];
		}
	}

	private int rand(int limit, int current)
	{
		int x;
		do
		{
			x = rand.nextInt(limit);
		}
		while (x == current);

		return x;
	}

	private void writeMatrix(String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write(header);
		for (int i = 0; i < matrix.length; i++)
		{
			writer.write("\n" + genes.get(i));
			for (int v : matrix[i])
			{
				writer.write("\t" + v);
			}
		}
		writer.close();
	}

	public void writeManyShufflings(String source, String outDir, int cnt) throws IOException
	{
		loadMatrix(source);
		for (int i = 0; i < cnt; i++)
		{
			shuffle();
			writeMatrix(outDir + "/DataMatrix-" + System.currentTimeMillis() + ".txt");
			Collections.shuffle(geneIndices);
		}
	}

	public static void main(String[] args) throws IOException
	{
		String source = "/home/babur/Documents/mutex/TCGA/PanCan/1/1/DataMatrix.txt";
		String typeSource = "/home/babur/Documents/mutex/TCGA/PanCan/SampleToType.txt";
		String outFile = "/home/babur/Documents/mutex/TCGA/PanCan-shuffled-10/1/1/DataMatrix.txt";
		System.out.println("outFile = " + outFile);

		AltMatrixShuffler shuffler = new AltMatrixShuffler();
		shuffler.loadMatrix(source);
		shuffler.loadSubtypes(typeSource);

		shuffler.shuffle();
		shuffler.writeMatrix(outFile);

//		shuffler.writeManyShufflings(source, "/home/babur/Documents/PanCan/ShuffledMatrices", 1000);
	}
}
