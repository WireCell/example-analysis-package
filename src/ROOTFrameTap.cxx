#include "WireCellExampleAna/ROOTFrameTap.h"
#include "WireCellIface/ITrace.h"
#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"

#include "WireCellAux/FrameTools.h"
#include "WireCellUtil/NamedFactory.h"

#include <string>
#include <vector>

/// macro to register name - concrete pair in the NamedFactory
/// @param NAME - used to configure node in JSON/Jsonnet
/// @parame CONCRETE - C++ concrete type
/// @parame ... - interfaces
WIRECELL_FACTORY(ROOTFrameTap, WireCell::Root::ROOTFrameTap,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;

Root::ROOTFrameTap::ROOTFrameTap() : log(Log::logger("ana")) {}

Root::ROOTFrameTap::~ROOTFrameTap() {}

void Root::ROOTFrameTap::configure(const WireCell::Configuration &cfg) {
  std::string fn;

  fn = cfg["output_filename"].asString();
  if (fn.empty()) {
    THROW(ValueError() << errmsg{
              "Must provide output filename to ROOTFrameTap"});
  }

  auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
  m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

  m_cfg = cfg;

  recreate_out_file();
}

WireCell::Configuration Root::ROOTFrameTap::default_configuration() const {
  Configuration cfg;

  cfg["anode"] = "AnodePlane";

  // Name of ROOT file to write.
  cfg["output_filename"] = "";

  cfg["nrebin"] = 1;

  return cfg;
}

namespace {
/// categorize traces to plane (u, v, w) via channel ID - plane ID relation
/// storaged in APA figure out time and channel binning for each plane
std::vector<WireCell::Binning> collate_byplane(const ITrace::vector &traces,
                                               const IAnodePlane::pointer anode,
                                               ITrace::vector byplane[]) {
  std::vector<int> uvwt[4];
  for (auto trace : traces) {
    const int chid = trace->channel();
    auto wpid = anode->resolve(chid);
    const int iplane = wpid.index();
    if (iplane < 0 || iplane >= 3) {
      THROW(RuntimeError() << errmsg{"Illegal wpid"});
    }
    uvwt[iplane].push_back(chid);
    byplane[iplane].push_back(trace);
    uvwt[3].push_back(trace->tbin());
    uvwt[3].push_back(trace->tbin() + trace->charge().size());
  }

  std::vector<Binning> binnings(4);
  for (int ind = 0; ind < 4; ++ind) {
    auto const &one = uvwt[ind];
    // std::cerr << "[wgu] get ind: " << ind << " size: " << one.size() <<
    // std::endl;
    if (one.empty()) {
      // THROW(ValueError() << errmsg{"ROOTFrameTap: bogus bounds"});
      std::cerr << "[wgu] plane: " << ind << " has not traces. " << std::endl;
    } else {
      auto mme = std::minmax_element(one.begin(), one.end());
      const int vmin = *mme.first;
      const int vmax = *mme.second;
      if (ind == 3) {
        const int n = vmax - vmin;
        // binnings.push_back(Binning(n, vmin, vmax););
        binnings.at(ind) = Binning(n, vmin, vmax);
      } else {
        // Channel-centered binning
        const double diff = vmax - vmin;
        // binnings.push_back(Binning(diff+1, vmin-0.5, vmax+0.5));
        binnings.at(ind) = Binning(diff + 1, vmin - 0.5, vmax + 0.5);
      }
    }
  }
  return binnings;
}
} // namespace

void Root::ROOTFrameTap::recreate_out_file() const {
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = "RECREATE";
  TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());
  output_tf->Close("R");
  delete output_tf;
  output_tf = nullptr;
}

void Root::ROOTFrameTap::peak_frame(const IFrame::pointer &frame) const {

  std::string frame_tags = "";
  for (auto tag : frame->frame_tags()) {
    frame_tags += " ";
    frame_tags += tag;
  }

  for (auto tag : frame->trace_tags()) {
    log->info(
        "ROOTFrameTap: trace tag: \"{}\" in frame: tags \"{}\" - ident: {}",
        tag, frame_tags, frame->ident());
  }
}

void Root::ROOTFrameTap::fill_hist(const IFrame::pointer &frame,
                                     TFile *output_tf) const {

  for (auto tag : frame->trace_tags()) {
    ITrace::vector traces = Aux::tagged_traces(frame, tag);
    if (traces.empty()) {
      log->warn("ROOTFrameTap: no tagged traces for \"{}\"", tag);
      continue;
    }

    log->debug("ROOTFrameTap: tag: \"{}\" with {} traces", tag,
               traces.size());

    ITrace::vector traces_byplane[3];
    auto binnings = collate_byplane(traces, m_anode, traces_byplane);
    Binning tbin = binnings[3];
    for (int iplane = 0; iplane < 3; ++iplane) {
      if (traces_byplane[iplane].empty())
        continue;
      const std::string name = Form("h%c_%s", 'u' + iplane, tag.c_str());
      Binning cbin = binnings[iplane];
      std::stringstream ss;
      ss << "ROOTFrameTap:"
         << " cbin:" << cbin.nbins() << "[" << cbin.min() << "," << cbin.max()
         << "]"
         << " tbin:" << tbin.nbins() << "[" << tbin.min() << "," << tbin.max()
         << "]";
      log->debug(ss.str());

      // consider to add nrebin ...
      int nbins = tbin.nbins();

      TH2F *hist =
          new TH2F(name.c_str(), name.c_str(), cbin.nbins(), cbin.min(),
                   cbin.max(), nbins, tbin.min(), tbin.max());

      hist->SetDirectory(output_tf);

      for (auto trace : traces_byplane[iplane]) {
        const int tbin1 = trace->tbin();
        const int ch = trace->channel();
        auto const &charges = trace->charge();
        for (size_t itick = 0; itick < charges.size(); ++itick) {

          int ibin = (tbin1 - tbin.min() + itick);

          hist->SetBinContent(
              cbin.bin(ch) + 1, ibin + 1,
              charges[itick] + hist->GetBinContent(cbin.bin(ch) + 1, ibin + 1));
        }
      }
    }
  }
}

IFrame::pointer
Root::ROOTFrameTap::fft_frame(const IFrame::pointer &in) const {
  WireCell::IFrame::trace_list_t indices;

  ITrace::vector *itraces = new ITrace::vector;

  for (auto trace : *(in->traces()).get()) {
    auto wf = trace->charge();
    auto wf_fft = WireCell::Waveform::magnitude(WireCell::Waveform::dft(wf));

    auto chid = trace->channel();
    auto tbin = trace->tbin();

    WireCell::SimpleTrace *trace_fft =
        new WireCell::SimpleTrace(chid, tbin, wf_fft);
    const size_t trace_index = itraces->size();
    indices.push_back(trace_index);
    itraces->push_back(ITrace::pointer(trace_fft));
  }

  SimpleFrame *sframe =
      new SimpleFrame(in->ident(), in->time(), ITrace::shared_vector(itraces),
                      in->tick(), in->masks());

  sframe->tag_frame("ana");

  sframe->tag_traces("ana_fft", indices);

  return IFrame::pointer(sframe);
}

bool Root::ROOTFrameTap::operator()(const IFrame::pointer &frame,
                                      IFrame::pointer &out) {

  out = frame;

  if (!frame) {
    // eos
    log->debug("ROOTFrameTap: EOS");
    return true;
  }
  if (frame->traces()->empty()) {
    log->debug("ROOTFrameTap: passing through empty frame ID {}",
               frame->ident());
    return true;
  }

  /// open output file
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = "UPDATE";
  log->debug("ROOTFrameTap: opening for output: {} with mode {}", ofname,
             mode);
  TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());

  peak_frame(frame);

  fill_hist(frame, output_tf);

  // auto frame_fft = fft_frame(frame);
  // fill_hist(frame_fft, output_tf);

  /// write and close output file
  auto count = output_tf->Write();
  log->debug("ROOTFrameTap: closing output file {}, wrote {} bytes", ofname,
             count);
  output_tf->Close();
  delete output_tf;
  output_tf = nullptr;
  return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
