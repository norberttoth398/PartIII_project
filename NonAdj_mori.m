function NonAdj_mori(grains, region, box, ideal_mori, threshold)

    %{
    The present function was created for non-adjacent grain misorientation
    analysis within a given region of interest using EBSD orientation
    maps. The function calculates misorentations between spatially
    constrained Ni3Al and NiAl grains, which is used to compare against an
    idealised misorientation given. This builds directly on MTEX 5.4.0, so 
    the environment needs to be initialised.

    INPUTS
    -------------------
    grains      - grains reconstructed from the enire map.
    region      - region for analysis to be perfomed on, this should be in
                vector format: [a b c d]
    box         - size of the box to be drawn around the Ni3Al grains for
                the spatial constraint.
    ideal_mori  - MTEX orientation object of the ideal misorientation
                considered here
    threshold   - Threshold used for the comparison between the
                misorientations calculated and the idealised one. (NUMBER,
                - int, float, double, etc. - not degrees type).

    OUTPUT
    -------------------
    A plot of the region of interesr considered here with grain pairs
    plotted colourcoded.
    %}
    
    %set up the grains to inspect
    grs = inpolygon(grains, region)
    grain = grains(grs)
    grains_region = grain(grain.grainSize > 2)

    %separate the two phases
    Ni3Al = grains_region('Ni3Al');
    NiAl = grains_region('NiAl');

    clear Relations
    for i = 1:length(Ni3Al) %iterate over all gamma' grains    

        for j = 1:length(NiAl) %iterate over all beta grains within 'box'
            in_box = inpolygon(NiAl(j),box);

            if in_box == 1
                mori_g = inv(NiAl(j).meanOrientation)*Ni3Al(i).meanOrientation;
                close = angle(ideal_mori, mori_g) < threshold*degree;

                if close == 1
                    Relations(i,j) = close_KS;
                else
                    Relations(i,j) = logical(0);
                end
            else
                Relations(i,j) = logical(0);
            end
        end
    end

    figure()
    plot(grains_region.boundary)
    hold on
    for i = 1:length(Ni3Al)
        condition = ismember(1, Relations(i,:));
        if condition == 1
            col = rand(1,3);
            plot(Ni3Al.subSet(i), 'micronBar', 'off', 'lineColor', col, 'linewidth', 1)
            hold on
            plot(NiAl.subSet(Relations(i,:)), 'micronbar', 'off', 'lineColor', col, 'linewidth', 1)
            hold on
        end
    end
    hold off
    legend(gca,'off')
    xlim([region(1) region(3)+region(1)])
    ylim([region(2) region(4)+region(2)])


    hold off
    legend(gca,'off')

end
