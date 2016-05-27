package org.panda.misc;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Pathway;
import org.panda.utility.ArrayUtil;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * For finding out the pathways related to the interactions in a SIF file, given that mediator IDs are provided in the
 * SIF file and the associated BioPAX file is available.
 *
 * @author Ozgun Babur
 */
public class MapSIFToPathways
{
	public static void main(String[] args) throws FileNotFoundException
	{
		printRelatedPathwayNames("/home/babur/Documents/Analyses/Sheep4x4/AR_vs_AC-paths-between.sif",
			"/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
	}

	public static void printRelatedPathwayNames(String siffile, String biopaxOwl) throws FileNotFoundException
	{
		SimpleIOHandler io = new SimpleIOHandler();
		Model model = io.convertFromOWL(new FileInputStream(biopaxOwl));
		Map<String, Set<String>> relToMeds = readMediators(siffile);
		Map<String, Set<String>> relToPN = mapMediatorsToPathwayNames(relToMeds, model);

		Map<String, Set<String>> patNameToRel = new HashMap<>();
		relToPN.keySet().forEach(rel ->
			relToPN.get(rel).forEach(name ->
			{
				if (!patNameToRel.containsKey(name)) patNameToRel.put(name, new HashSet<>());
				patNameToRel.get(name).add(rel);
			}));

		patNameToRel.keySet().stream()
			.sorted((p1, p2) -> new Integer(patNameToRel.get(p2).size()).compareTo(patNameToRel.get(p1).size()))
			.forEach(p ->
			{
				System.out.print(p + "\t");
				for (String rel : patNameToRel.get(p))
				{
					System.out.print(rel + ", ");
				}
				System.out.println();
			});
	}

	private static Map<String, Set<String>> readMediators(String siffile) { try
	{
		Map<String, Set<String>> map = new HashMap<>();

		Files.lines(Paths.get(siffile)).map(l -> l.split("\t")).filter(t -> t.length > 3).forEach(token ->
			map.put(ArrayUtil.getString(" ", token[0], token[1], token[2]),
				new HashSet<>(Arrays.asList(token[3].split(";")))));

		return map;
	}
	catch (Exception e){throw new RuntimeException(e);}}

	private static Set<String> getPathwayNames(Set<BioPAXElement> mediators)
	{
		PathAccessor pa1 = new PathAccessor("Interaction/stepProcessOf/pathwayOrderOf");
		PathAccessor pa2 = new PathAccessor("Interaction/pathwayComponentOf");

		Set<Pathway> pathways = new HashSet<>();
		pathways.addAll(pa1.getValueFromBeans(mediators));
		pathways.addAll(pa2.getValueFromBeans(mediators));

		return pathways.stream().map(Pathway::getDisplayName).collect(Collectors.toSet());
	}

	private static Set<BioPAXElement> findBioPAXElements(Set<String> mediators, Model model)
	{
		return mediators.stream().map(model::getByID).collect(Collectors.toSet());
	}

	private static Map<String, Set<String>> mapMediatorsToPathwayNames(Map<String, Set<String>> relToMeds, Model model)
	{
		Map<String, Set<String>> relToPN = new HashMap<>();

		relToMeds.keySet().forEach(rel ->
		{
			Set<BioPAXElement> eles = findBioPAXElements(relToMeds.get(rel), model);
			Set<String> names = getPathwayNames(eles);
			if (!names.isEmpty()) relToPN.put(rel, names);
		});

		return relToPN;
	}
}
