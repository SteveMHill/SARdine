// tests/annotation_parser.rs

#![allow(clippy::approx_constant)]

fn approx(a: f64, b: f64, eps: f64) {
    assert!((a - b).abs() <= eps, "approx failed: |{a} - {b}| > {eps}");
}

mod s1 {
    use super::approx;

    use crate as _ignore; // keep cargo happy when running as workspace
    use std::env;

    // Pull the parser + helpers from your crate
    use sardine::io::annotation::{parse_annotation_xml, AnnotationRoot as ProductRoot};

    // Handy constant
    const C: f64 = sardine::constants::physical::SPEED_OF_LIGHT_M_S as f64;

    // ---------- Fixtures ----------

    // Full-ish SAFE-like annotation with namespaces stripped by the parser,
    // includes:
    // - generalAnnotation/productInformation (pixel spacing, radar freq, etc)
    // - downlinkInformationList/prf (so PRF beats azimuthFrequency)
    // - imageAnnotation/imageInformation (dims, azimuthFrequency)
    // - swathProcParamsList for IW1/2/3 (subswaths present)
    // - swathTiming with burst, valid sample arrays, large byteOffset (u64)
    // - geolocationGrid (for bbox + near/far)
    // - antennaPattern with space-separated arrays
    // - SAFE-style <orbit> blocks with nested position/velocity
    const XML_SAFE: &str = r#"
<?xml version="1.0" encoding="UTF-8"?>
<s1:product xmlns:s1="s1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <adsHeader>
    <missionId>S1A</missionId>
    <productType>SLC</productType>
    <startTime>2020-12-28T21:59:42.123456Z</startTime>
    <stopTime>2020-12-28T21:59:52Z</stopTime>
    <absoluteOrbitNumber>30639</absoluteOrbitNumber>
    <missionDataTakeId>123456</missionDataTakeId>
    <imageNumber>1</imageNumber>
  </adsHeader>

  <generalAnnotation>
    <productInformation>
      <platformHeading>12.34</platformHeading>
      <rangeSamplingRate>4.991e+07</rangeSamplingRate>
      <radarFrequency>5.405e+09</radarFrequency>
      <azimuthSteeringRate>0.001</azimuthSteeringRate>
      <rangePixelSpacing>2.329560</rangePixelSpacing>
      <azimuthPixelSpacing>13.943035</azimuthPixelSpacing>
    </productInformation>

    <downlinkInformationList>
      <downlinkInformation>
        <prf>1710.0</prf>
      </downlinkInformation>
    </downlinkInformationList>

    <dcEstimateList>
      <dcEstimate>
        <t0>0.0</t0>
        <dataDcPolynomial>100.0 2.0 -0.5</dataDcPolynomial>
      </dcEstimate>
    </dcEstimateList>
  </generalAnnotation>

  <imageAnnotation>
    <imageInformation>
      <slantRangeTime>0.004</slantRangeTime>
      <rangePixelSpacing>2.33</rangePixelSpacing>
      <azimuthPixelSpacing>13.9</azimuthPixelSpacing>
      <numberOfSamples>100</numberOfSamples>
      <numberOfLines>200</numberOfLines>
      <azimuthTimeInterval>0.00028</azimuthTimeInterval>
      <azimuthFrequency>1725.0</azimuthFrequency>
      <productFirstLineUtcTime>2020-12-28T21:59:42.123456Z</productFirstLineUtcTime>
      <productLastLineUtcTime>2020-12-28T21:59:52Z</productLastLineUtcTime>
    </imageInformation>

    <processingInformation>
      <swathProcParamsList>
        <swathProcParams><swath>IW1</swath></swathProcParams>
        <swathProcParams><swath>IW2</swath></swathProcParams>
        <swathProcParams><swath>IW3</swath></swathProcParams>
      </swathProcParamsList>
    </processingInformation>
  </imageAnnotation>

  <swathTiming>
    <burstList>
      <burst>
        <azimuthTime>2020-12-28T21:59:42Z</azimuthTime>
        <azimuthAnxTime>1234.5</azimuthAnxTime>
        <sensingTime>2020-12-28T21:59:42Z</sensingTime>
        <byteOffset>5000000000</byteOffset>
        <firstValidSample>10 10 10</firstValidSample>
        <lastValidSample>90 90 90</lastValidSample>
      </burst>
    </burstList>
  </swathTiming>

  <geolocationGrid>
    <geolocationGridPointList>
      <geolocationGridPoint>
        <slantRangeTime>0.003</slantRangeTime>
        <line>0</line><pixel>0</pixel>
        <latitude>10.0</latitude><longitude>20.0</longitude><height>0</height>
        <incidenceAngle>30.0</incidenceAngle><elevationAngle>5.0</elevationAngle>
      </geolocationGridPoint>
      <geolocationGridPoint>
        <slantRangeTime>0.005</slantRangeTime>
        <line>1</line><pixel>1</pixel>
        <latitude>12.0</latitude><longitude>22.0</longitude><height>0</height>
        <incidenceAngle>40.0</incidenceAngle><elevationAngle>6.0</elevationAngle>
      </geolocationGridPoint>
    </geolocationGridPointList>
  </geolocationGrid>

  <antennaPattern>
    <antennaPattern>
      <elevationAngle>1 2 3</elevationAngle>
      <incidenceAngle>30 35 40</incidenceAngle>
      <slantRangeTime>0.003 0.004 0.005</slantRangeTime>
      <elevationPattern>1 2 3</elevationPattern>
      <terrainHeight>0 0 0</terrainHeight>
      <roll>0</roll>
    </antennaPattern>
  </antennaPattern>

  <!-- SAFE-style orbit blocks -->
  <orbit>
    <time>2020-12-28T21:59:42.123456Z</time>
    <position><x>1</x><y>2</y><z>3</z></position>
    <velocity><x>7000</x><y>0</y><z>0</z></velocity>
  </orbit>
  <orbit>
    <time>2020-12-28T21:59:43Z</time>
    <position><x>2</x><y>3</y><z>4</z></position>
    <velocity><x>7000</x><y>10</y><z>0</z></velocity>
  </orbit>
</s1:product>
"#;

    // ZIP-like fixture:
    // - No downlinkInformationList/prf (forces PRF fallback to azimuthFrequency)
    // - ZIP-style <stateVector> block with vx/vy/vz
    const XML_ZIP_STYLE: &str = r#"
<product>
  <adsHeader>
    <missionId>S1B</missionId>
    <startTime>2020-12-28T21:59:42+00:00</startTime>
    <stopTime>2020-12-28T21:59:52+00:00</stopTime>
    <absoluteOrbitNumber>30640</absoluteOrbitNumber>
    <missionDataTakeId>123457</missionDataTakeId>
    <imageNumber>1</imageNumber>
  </adsHeader>

  <generalAnnotation>
    <productInformation>
      <rangeSamplingRate>4.991e+07</rangeSamplingRate>
      <radarFrequency>5.405e+09</radarFrequency>
      <rangePixelSpacing>2.5</rangePixelSpacing>
      <azimuthPixelSpacing>14.0</azimuthPixelSpacing>
    </productInformation>
  </generalAnnotation>

  <imageAnnotation>
    <imageInformation>
      <azimuthFrequency>1700.0</azimuthFrequency>
      <slantRangeTime>0.0045</slantRangeTime>
      <numberOfSamples>50</numberOfSamples>
      <numberOfLines>60</numberOfLines>
      <productFirstLineUtcTime>2020-12-28T21:59:42Z</productFirstLineUtcTime>
      <productLastLineUtcTime>2020-12-28T21:59:52Z</productLastLineUtcTime>
    </imageInformation>
    <processingInformation>
      <swathProcParamsList>
        <swathProcParams><swath>IW2</swath></swathProcParams>
      </swathProcParamsList>
    </processingInformation>
  </imageAnnotation>

  <geolocationGrid>
    <geolocationGridPointList>
      <geolocationGridPoint>
        <slantRangeTime>0.0041</slantRangeTime>
        <line>0</line><pixel>0</pixel>
        <latitude>0</latitude><longitude>0</longitude><height>0</height>
        <incidenceAngle>35</incidenceAngle><elevationAngle>0</elevationAngle>
      </geolocationGridPoint>
      <geolocationGridPoint>
        <slantRangeTime>0.0049</slantRangeTime>
        <line>1</line><pixel>1</pixel>
        <latitude>1</latitude><longitude>1</longitude><height>0</height>
        <incidenceAngle>45</incidenceAngle><elevationAngle>0</elevationAngle>
      </geolocationGridPoint>
    </geolocationGridPointList>
  </geolocationGrid>

  <stateVector>
    <time>2020-12-28T21:59:42Z</time>
    <x>7000000</x><y>0</y><z>0</z>
    <vx>0</vx><vy>7600</vy><vz>0</vz>
  </stateVector>
</product>
"#;

    // ---------- Tests ----------

    #[test]
    fn parse_basic_and_core_fields_safe() {
        let parsed = parse_annotation_xml(XML_SAFE).expect("parse safe");
        // Pixel spacing (prefer generalAnnotation)
        let (rps, aps) = parsed.get_pixel_spacing().expect("pixel spacing");
        approx(rps, 2.329560, 1e-9);
        approx(aps, 13.943035, 1e-6);

        // PRF should come from downlink list (not azimuthFrequency)
        let prf = parsed.get_pulse_repetition_frequency().expect("prf");
        approx(prf, 1710.0, 1e-9);

        // Radar frequency + wavelength
        let rf = parsed.get_radar_frequency_hz().expect("radar freq");
        approx(rf, 5.405e9, 1e3);
        let lambda = (sardine::constants::physical::SPEED_OF_LIGHT_M_S as f64) / rf;
        approx(lambda, 0.0555, 5e-4);

        // Slant range time from imageInformation
        let srt = parsed.get_slant_range_time().expect("srt");
        approx(srt, 0.004, 1e-12);
    }

    #[test]
    fn bounding_box_from_geogrid() {
        let parsed = parse_annotation_xml(XML_SAFE).expect("parse");
        let bbox = ProductRoot::extract_bounding_box(&parsed).expect("bbox");
        approx(bbox.min_lat, 10.0, 1e-9);
        approx(bbox.max_lat, 12.0, 1e-9);
        approx(bbox.min_lon, 20.0, 1e-9);
        approx(bbox.max_lon, 22.0, 1e-9);
    }

    #[test]
    fn antenna_pattern_space_separated_arrays() {
    let parsed = parse_annotation_xml(XML_SAFE).expect("parse");
    let sets = parsed
      .get_antenna_patterns()
      .expect("antenna pattern sets");
    assert_eq!(sets.len(), 1);
    let (elev, inc, srt) = &sets[0];
        assert_eq!(elev, &vec![1.0, 2.0, 3.0]);
        assert_eq!(inc, &vec![30.0, 35.0, 40.0]);
        assert_eq!(srt, &vec![0.003, 0.004, 0.005]);
    }

    #[test]
    fn burst_arrays_and_byteoffset_u64() {
        let parsed = parse_annotation_xml(XML_SAFE).expect("parse");
        let st = parsed.swath_timing.as_ref().expect("swathTiming");
        let bl = st.burst_list.as_ref().expect("burstList");
        let b = bl.bursts.as_ref().expect("bursts");
        assert_eq!(b.len(), 1);
        let burst = &b[0];

        // arrays parsed
        assert_eq!(burst.first_valid_sample, vec![10, 10, 10]);
        assert_eq!(burst.last_valid_sample, vec![90, 90, 90]);

        // large offset survives as u64 (for SAR files > 4GB)
        assert_eq!(burst.byte_offset, Some(5_000_000_000u64));
    }

    #[test]
    fn orbit_extraction_safe_style() {
        // Allow regex orbit enrichment (default true in suggested code).
        env::set_var("SARDINE_ALLOW_REGEX_ORBITS", "1");
        let parsed = parse_annotation_xml(XML_SAFE).expect("parse");
        let orbits = parsed.orbit_list.as_ref().expect("orbit list");
        assert_eq!(orbits.len(), 2);
        // velocity magnitudes ~7000..7010
        for sv in orbits {
            let v =
                (sv.velocity[0].powi(2) + sv.velocity[1].powi(2) + sv.velocity[2].powi(2)).sqrt();
            assert!(v > 6900.0 && v < 7050.0, "v={v}");
        }
    }

    #[test]
    fn orbit_extraction_zip_style_and_prf_fallback() {
        env::set_var("SARDINE_ALLOW_REGEX_ORBITS", "1");
        let parsed = parse_annotation_xml(XML_ZIP_STYLE).expect("parse");
        // PRF should fallback to azimuthFrequency when downlink list is absent
        let prf = parsed
            .get_pulse_repetition_frequency()
            .expect("prf fallback");
        approx(prf, 1700.0, 1e-9);

        // ZIP-style <stateVector> present
        let orbits = parsed.orbit_list.as_ref().expect("orbit list");
        assert_eq!(orbits.len(), 1);
        approx(orbits[0].position[0], 7_000_000.0, 1.0);
        approx(orbits[0].velocity[1], 7_600.0, 0.1);
    }

    #[test]
    fn doppler_centroid_polynomial_evaluation() {
        let parsed = parse_annotation_xml(XML_SAFE).expect("parse");
        // t0 = 0, coeffs = [100, 2, -0.5]: f(1) = 100 + 2*1 + -0.5*1^2 = 101.5
        let f_at_1 = parsed.evaluate_doppler_centroid(1.0).expect("doppler");
        approx(f_at_1, 101.5, 1e-9);
    }

    #[test]
    fn range_doppler_params_are_consistent() {
        let parsed = parse_annotation_xml(XML_SAFE).expect("parse");
        let params = parsed.extract_range_doppler_params().expect("rd params");

        approx(params.range_pixel_spacing, 2.329560, 1e-9); // from generalAnnotation
        approx(params.azimuth_pixel_spacing, 13.943035, 1e-6);
        approx(params.slant_range_time, 0.004, 1e-12); // from imageInformation
        approx(params.prf, 1710.0, 1e-9);
        approx(params.wavelength, C / 5.405e9, 1e-6);
        assert!(params.product_start_time_abs > 1_600_000_000.0);
        // doppler model present
        let dc = params.doppler_centroid.as_ref().expect("dc model");
        assert_eq!(dc.coeffs.len(), 3);
    }

    #[test]
    fn subswaths_present_for_all_and_specific() {
        let parsed_safe = parse_annotation_xml(XML_SAFE).expect("parse");
        let map = ProductRoot::extract_subswaths(&parsed_safe).expect("subs");
        assert!(map.contains_key("IW1"));
        assert!(map.contains_key("IW2"));
        assert!(map.contains_key("IW3"));

        let parsed_zip = parse_annotation_xml(XML_ZIP_STYLE).expect("parse");
        let map2 = ProductRoot::extract_subswaths(&parsed_zip).expect("subs2");
        assert!(map2.contains_key("IW2"));
        assert_eq!(map2.len(), 1);
    }

    #[test]
    fn time_parsing_variants_survive_via_rd_params() {
        // Same doc but vary start/stop time formats; reuse XML_ZIP_STYLE as template and patch the times.
        let variants = [
            ("2020-12-28T21:59:42+00:00", "2020-12-28T21:59:52+00:00"),
            (
                "2020-12-28T21:59:42.123456+00:00",
                "2020-12-28T21:59:52.5+00:00",
            ),
            ("2020-12-28T21:59:42Z", "2020-12-28T21:59:52Z"),
            ("2020-12-28T21:59:42.123456Z", "2020-12-28T21:59:52.999999Z"),
        ];

        for (start, stop) in variants {
            let xml = XML_ZIP_STYLE
                .replace("2020-12-28T21:59:42+00:00", start)
                .replace("2020-12-28T21:59:52+00:00", stop);
            let parsed = parse_annotation_xml(&xml).expect("parse");
            // Will parse product start time internally; we just ensure it does not error.
            let params = parsed.extract_range_doppler_params().expect("rd params");
            assert!(params.product_start_time_abs > 1.0e9);
        }
    }
}
