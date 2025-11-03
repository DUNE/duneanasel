#pragma once

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

template <typename SRTypeConstraint, typename T> struct Proxyable {

  constexpr static bool is_type = std::is_same_v<SRTypeConstraint, T>;
  constexpr static bool is_proxy_to_type =
      std::is_same_v<typename caf::Proxy<SRTypeConstraint>, T>;

  static_assert(is_type || is_proxy_to_type,
                "Invalid type used for template function");

  using type = std::enable_if_t<is_type || is_proxy_to_type, T>;
};

template <typename SRType, typename T>
using Proxyable_t = typename Proxyable<SRType, T>::type;
