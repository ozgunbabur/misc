package org.panda.misc;

import org.panda.utility.CollectionUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by babur on 4/22/16.
 */
public class MutexReader
{
	public static Set<Group> readMutexResults(String dir)
	{try{
		String file = dir + "/ranked-groups.txt";
		boolean hasQval;
		hasQval = Files.lines(Paths.get(file)).limit(1).anyMatch(line -> line.split("\t").length >= 3);

		Set<Group> groups = new HashSet<>();
		Files.lines(Paths.get(file)).skip(1).map(line -> line.split("\t"))
			.filter(token -> Double.parseDouble(token[0]) < 1).forEach(token ->
				groups.add(new Group(Arrays.asList(token).subList(hasQval ? 2 : 1, token.length),
					dir, Double.parseDouble(token[0]))));

		return groups;
	} catch (IOException e){throw new RuntimeException(e);}}

	public static Set<Group> readMutexResultsRecursive(String dir, Set<Group> result)
	{
		if (new File(dir + "/ranked-groups.txt").exists())
			result.addAll(readMutexResults(dir));

		Arrays.stream(new File(dir).listFiles()).filter(File::isDirectory).forEach(d ->
			readMutexResultsRecursive(d.getPath(), result));

		return result;
	}

	static class Group
	{
		List<String> genes;
		String fromDir;
		double score;

		public Group(List<String> genes, String fromDir, double score)
		{
			this.genes = genes;
			this.fromDir = fromDir;
			this.score = score;
		}

		public double getScore()
		{
			return score;
		}

		@Override
		public String toString()
		{
			return score + "\t" + CollectionUtil.merge(genes, " ") + "\t" + fromDir;
		}

		public void shortenLoc(String remove)
		{
			if (fromDir.startsWith(remove)) fromDir = fromDir.substring(remove.length());
		}

		@Override
		public int hashCode()
		{
			int h = 0;
			for (String gene : genes)
			{
				h += gene.hashCode();
			}
			h += fromDir.hashCode();
			h += Double.hashCode(score);
			return h;
		}

		@Override
		public boolean equals(Object obj)
		{
			return obj instanceof Group && ((Group) obj).genes.containsAll(genes) &&
				((Group) obj).genes.size() == genes.size() && ((Group) obj).fromDir.equals(fromDir) &&
				((Group) obj).score == score;
		}
	}

	public static void main(String[] args)
	{
		String dir = "/home/babur/Documents/mutex/TCGA/";
		readMutexResultsRecursive(dir, new HashSet<>()).stream()
			.sorted(Comparator.comparing(Group::getScore)).filter(g -> g.score < 0.05).peek(g -> g.shortenLoc(dir))
			.map(g -> g.toString().replaceAll("/", "\t")).forEach(System.out::println);
	}
}
